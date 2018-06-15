/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* Last Modified: 12APR2009 */

/* De-macro'd by the ACE-team 18MAY2012 wattup! */

/*Converted into template classes by CTS 11DEC2014*/

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <cctype>
#include <string>
#include <cstdlib>
#include <unistd.h>

#if HAVE_ZLIB
#include <zlib.h>
class FunctorZlib 
{
public:
    int operator()(gzFile file, void * buffer, unsigned int len)
    {
        return gzread(file, buffer, len);
    }
};
#endif


#if HAVE_BZIP2
#include <bzlib.h>
class FunctorBZlib2
{
public:
    int operator()(BZFILE* file, void * buffer, int len) 
    {
        return BZ2_bzread(file, buffer, len );
    }
};
#endif

class FunctorRead
{
public:
    ssize_t operator()(int fd, void *buf, size_t count)
    {
        return read(fd, buf, count);
    }
};

class kseq
{
public:
    ~kseq() = default;
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;
    int last_char = 0;
};


template<class ret_t, class ReadFunction>
class kstream
{
public:
    enum Delimiter {SPACE=1, EOL=0};
    kstream(ret_t f, ReadFunction rf) : bufferSize(4096), f(f), is_eof(0), begin(0), end(0), readfunc(rf) {
        buf = (char*) malloc(bufferSize);
        c = getHeader();
    }

    ~kstream()
    {
        free(buf);
    }

    int readFastq(kseq& seq) {
        seq.comment.clear();
        seq.seq.clear();
        seq.qual.clear();

        if (c != '\n') {
            std::cerr << "ERROR: Unexpected record start, last valid record: " << seq.name << std::endl;
            return -3;
        }

        if (getName(seq) == -1) {
            if (is_eof != 1) {
                std::cerr << "ERROR: Unexpected ID after " << seq.name << std::endl;
                return -3;
            } else {
                return -2;
            }
        }

        bool good(false);
        switch(c){
            case 10:    // new line character
                good = true;
                break;
        }
        if (!good) {
            std::cerr << "ERROR: There must have been a problem reading the ID of " << seq.name << std::endl;
            return -3;
        }

        getSeq(seq);
        if (c != '+')
            return (int)seq.seq.length();
        while ((c = this->getc()) != -1 && c != '\n');  // Ignore whatever comes after '+'

        if (c == -1) {
            std::cerr << "ERROR: File ended unexpectedly on ID " << seq.name << std::endl;
            return -3; // File ended unexpectedly
        }

        getQual(seq);
        if (c != '\n') {
            std::cerr << "ERROR: The quality string appears to be longer for ID " << seq.name << std::endl;
            return -3;
        }
        if (seq.seq.length() != seq.qual.length()) {
            std::cerr << "ERROR: The length of the sequence doesn't match the quality for ID " << seq.name << std::endl;
            return -3;
        }
        return (int)seq.seq.length();
    }

    int readFasta(kseq& seq) {
        seq.comment.clear();
        seq.seq.clear();
        seq.qual.clear();

        if (c != '\n') {
            std::cerr << "ERROR: Unexpected record start, last valid record: " << seq.name << std::endl;
            return -3;
        }

        if (getName(seq) == -1) {
            if (is_eof != 1) {
                std::cerr << "ERROR: Unexpected ID after " << seq.name << std::endl;
                return -3;
            } else {
                return -2;
            }
        }

        bool good(false);
        switch (c) {
            case 10:    // new line character
                good = true;
                break;
        }
        if (!good) {
            std::cerr << "ERROR: There must have been a problem reading the ID of " << seq.name << std::endl;
            return -3;
        }

        getSeq(seq);
        if (c == '\n') {
            return (int) seq.seq.length();
        } else if (is_eof != 1) {
            std::cerr << "ERROR: There was an error reading the sequence of the record: " << seq.name << std::endl;
            return -3;
        }
    }

private:
    int getHeader() {
        int c;
        while ((c = this->getc()) != -1 && c != '>' && c != '@');
        if (c == -1) {
            std::cerr << "ERROR: The file appears to be empty" << std::endl;
            return -1;
        }
        if (c == int('>') || c == int('@')){
            begin--;
            c = int('\n');
        }
        return c;
    }

    int getName(kseq& seq) {
        if (this->getuntil(SPACE, seq.name, &c) < 0)
            return -1;
        if (c != '\n')
            return this->getuntil( EOL, seq.comment, &c);
    }

    int getSeq(kseq& seq){
        while ((c = this->getc()) != -1 && c != '>' && c != '+' && c != '@')
        {
            if (isgraph(c))
            {
                seq.seq += (char)c;
            }
        }
    }

    void getQual(kseq& seq){
        while ((c = this->getc()) != -1 && seq.qual.length() < seq.seq.length()) {
            if (c >= 33 && c <= 127)
                seq.qual += (char)c;
        }
    }

    int getc()
    {
        if (this->is_eof && this->begin >= this->end)
            return -1;
        if (this->begin >= this->end)
        {
            this->begin = 0;
            this->end = this->readfunc(this->f, this->buf, bufferSize);
            if (this->end < bufferSize)
                this->is_eof = 1;
            if (this->end == 0)
                return -1;
        }
        return (int)this->buf[this->begin++];
    }

    int getuntil(Delimiter delimiter, std::string &str, int *dret)
    {
        if (dret)
            *dret = 0;
        if (!str.empty()) {
            str.clear();
        }

        if (this->begin >= this->end && this->is_eof)
            return -1;
        for (;;)
        {
            int i;
            if (this->begin >= this->end)
            {
                if (!this->is_eof)
                {
                    this->begin = 0;
                    this->end = this->readfunc(this->f, this->buf, bufferSize);
                    if (this->end < bufferSize)
                        this->is_eof = 1;
                    if (this->end == 0)
                        break;
                }
                else
                    break;
            }
            if (delimiter == EOL)
            {
                for (i = this->begin; i < this->end; ++i)
                {
                    if (this->buf[i] == '\n')
                        break;
                }
            }
            else if (delimiter == SPACE)
            {
                for (i = this->begin; i < this->end; ++i)
                {
                    if (isspace(this->buf[i]))
                        break;
                }
            }
            else i = 0;

            str.append(this->buf + this->begin, static_cast<unsigned long>(i - this->begin));
            this->begin = i + 1;
            if (i < this->end)
            {
                if (dret)
                    *dret = this->buf[i];
                break;
            }
        }
        return (int)str.length();
    }

    int c;
    char *buf;
    int begin;
    int end;
    int is_eof;
    ret_t f;
    ReadFunction readfunc;
    const unsigned int bufferSize;
};

#endif

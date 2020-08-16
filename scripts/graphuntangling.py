def solve_bubble(t,s,ge):
    #This checks if a bubble can be due to a sequencing error and queues a node deletion on the error node
    fapath=s.stride_out_in_order(t.frontiers[0].rc().node_id()).nodes[:3]
    fbpath=[-x for x in s.stride_out_in_order(t.frontiers[1].rc().node_id()).nodes[:3][::-1]]
    #print(fapath,fbpath)
    if len(fapath)<2: fapath=fbpath
    if len(fbpath)<2: fbpath=fapath
    if len(fapath)<2: return False
    if fapath==fbpath and abs(t.internals[0].node_id())==abs(fapath[1]):
        return ge.queue_node_deletion(t.internals[1].node_id())
    if fapath==fbpath and abs(t.internals[1].node_id())==abs(fapath[1]):
        return ge.queue_node_deletion(t.internals[0].node_id())

def solve_repeat(t,s,ge):
    #solves a single-node repeat
    rnid=t.internals[0].node_id()
    ins=[x.node() for x in t.internals[0].prev()]
    inids=[x.node_id() for x in ins]
    outs=[x.node() for x in t.internals[0].next()]
    onids=[x.node_id() for x in outs]
    solutions=[]
    if len(ins)!=len(outs): return False
    for inv in ins:
        fpath=s.stride_out_in_order(inv.node_id()).nodes[:3]
        if len(fpath)==3 and fpath[1]==rnid and fpath[2] in onids:
            bpath=[-x for x in s.stride_out_in_order(-fpath[2]).nodes[:3][::-1]]
            if bpath==fpath: solutions.append(fpath)
    if len(solutions)!=len(ins): return False
    ge.queue_node_expansion(rnid,[[-x[0],x[2]] for x in solutions])
    return True

def solve_tip(t,s,ge):
    #just check that each frontier connects to the other one through reads
    f0s=s.stride_out_in_order(-t.frontiers[0].node_id()).nodes
    f1s=s.stride_out_in_order(-t.frontiers[1].node_id()).nodes
    if len(f0s)>1 and f0s[1]==t.frontiers[1].node_id() and len(f1s)>1 and f1s[1]==t.frontiers[0].node_id():
        ge.queue_node_deletion(t.internals[0].node_id())
        return True
    return False

def end_to_end_solution(p,fnids):
    for i in range(len(p)):
        if -p[i] in fnids:
            for j in range(i+1,len(p)):
                if p[j] in fnids:
                    if abs(p[i])<abs(p[j]): sol=p[i:j+1]
                    else: sol=[-x for x in p[i:j+1][::-1]]
                    #go through sol, if there's a single-option missing node, add it
                    csol=[sol[0]]
                    for n in sol[1:]:
                        if n in ws.sdg.fw_reached_nodes(csol[-1],1):
                            csol.append(n)
                        elif n in ws.sdg.fw_reached_nodes(csol[-1],2):
                            a=[]
                            for skipped in ws.sdg.fw_reached_nodes(csol[-1],1):
                                if n in ws.sdg.fw_reached_nodes(skipped,1):
                                    a.append(skipped)
                            if len(a)==1:
                                csol.append(a[0])
                                csol.append(n)
                            else:
                                return sol
                        else:
                            return sol
                    return csol
            return []
    return []

def solve_unclassified(t,s,ge):
    #first stride from all frontiers:
    fs={}
    fnids=[x.node_id() for x in t.frontiers]
    sols={x:[] for x in fnids}
    for f in t.frontiers:
        fs[-f.node_id()]=s.stride_out_in_order(-f.node_id()).nodes
        #for x in fs: print(x,fs[x])
    ns={}
    for nv in t.internals:
        ns[nv.node_id()]=[-x for x in s.stride_out_in_order(-nv.node_id()).nodes[:1:-1]]+s.stride_out_in_order(nv.node_id()).nodes
    #for x in ns: print(x,ns[x])
    sols=[]
    for x in list(fs.values())+list(ns.values()):
        sol=end_to_end_solution(x,fnids)
        if sol and sol not in sols:
            sols.append(sol)
    #for x in sols: print(x)
    if len(sols)==len(fnids)/2 and len(set([x[0] for x in sols]+[-x[-1] for x in sols]))==len(fnids):
        #print("UNCLASSIFIED SOLVED:")
        for x in sols:
            try:
                ps=SDG.SequenceDistanceGraphPath(ws.sdg,x).sequence()
            except:
                return False
        for x in sols: ge.queue_path_detachment(x,True)
        for x in t.internals: ge.queue_node_deletion(x.node_id())
        return True
    return False

def solve_all_tangles(ws,fsize=220,fminkci=-1,fmaxkci=-1,apply=False):
    total=solved=unsolved=0
    s=SDG.Strider(ws)
    s.add_datastore(peds)
    ge=SDG.GraphEditor(ws)
    solved=Counter()
    unsolved=Counter()
    for t in ws.sdg.get_all_tangleviews(fsize,fminkci,fmaxkci):
        tc=classify_tangle(t)
        if tc=="tip":
            #TODO add tip to reconnection list
            if solve_tip(t,s,ge): solved[tc]+=1
            else: unsolved[tc]+=1
        elif tc=="bubble":
            if solve_bubble(t,s,ge): solved[tc]+=1
            else: unsolved[tc]+=1
        elif tc[:6]=="repeat":
            if solve_repeat(t,s,ge): solved[tc]+=1
            else: unsolved[tc]+=1
        else:
            if solve_unclassified(t,s,ge): solved[tc]+=1
            else: unsolved[tc]+=1
        total+=1
    print("Total tangles: %d " % (total) )
    print("Solved:   ",solved.most_common())
    print("Unsolved: ",unsolved.most_common())
    if apply:
        ge.apply_all()
        ws.sdg.join_all_unitigs()



def classify_tangle(t):
    if len(t.frontiers)==0:
        return "debris"
    if len(t.frontiers)==1:
        nids=set([x.node_id() for x in t.internals])
        seen_nids=set()
        nexts=[t.frontiers[0].rc()]
        abort=False
        while nexts and not abort:
            new_nexts=[]
            for nv in nexts:
                for nnl in nv.next():
                    if abs(nnl.node().node_id()) not in nids:
                        abort=True
                        break
                    if nnl.node().node_id() not in seen_nids:
                        seen_nids.add(nnl.node().node_id())
                        new_nexts.append(nnl.node())
                if abort: break
            nexts=new_nexts
        if not abort: return "complex tip"
    if len(t.internals)==1:
        nv=t.internals[0]
        if len(t.frontiers)==2 and len(nv.prev())+len(nv.next())==1:
            return "tip"
        if len(t.frontiers)==2*len(nv.prev()) and len(nv.prev())==len(nv.next()) and len(set([abs(nv.node_id())]+[abs(x.node().node_id()) for x in nv.prev()+nv.next()]))==len(nv.prev())*2+1:
            return "repeat %d:%d"%(len(nv.prev()),len(nv.prev()))
    if len(t.internals)==2 and len(t.frontiers)==2:
        i0p=t.internals[0].parallels()
        if len(i0p)==1 and i0p[0]==t.internals[1]:
            return "bubble"
        else:
            i0n=t.internals[0].next()
            i0p=t.internals[0].prev()
            i1n=t.internals[1].next()
            i1p=t.internals[1].prev()
            if len(i0n)==2 and len(i0p)==2 and len(i1n)==1 and len(i1p)==1 and i1n[0].node()==i1p[0].node():
                return "loop"
            if len(i1n)==2 and len(i1p)==2 and len(i0n)==1 and len(i0p)==1 and i0n[0].node()==i0p[0].node():
                return "loop"
    return "unclassified"
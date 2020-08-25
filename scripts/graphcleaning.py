import SDGpython as SDG
from collections import Counter

def is_canonical_repeat(nv):
    if len(nv.prev())==len(nv.next())>1:
        for pl in nv.prev():
            if len(pl.node().next())>1: return False
        for nl in nv.next():
            if len(nl.node().prev())>1: return False
        return True
    return False

def is_bubble_side(nv):
    if len(nv.parallels())==len(nv.prev())==len(nv.next())==1 and len(nv.prev()[0].node().next())==len(nv.next()[0].node().prev())==2: return True
    return False

def is_tip(nv):
    if (len(nv.prev())==0 and len(nv.next())>0) or (len(nv.next())==0 and len(nv.prev())>0): return True
    return False

def simple_structures_stats(dg):
    total=canrep=bubsides=tips=0

    for nv in dg.get_all_nodeviews():
        total+=1
        if is_canonical_repeat(nv): canrep+=1
        if is_bubble_side(nv): bubsides+=1
        if is_tip(nv): tips+=1
    print("%d tips, %d canonical repeats and %d bubble sides out of %d nodes"%(tips,canrep,bubsides,total))


def solve_canonical_repeat(ge,nv,peds,min_support=5,max_noise=10,snr=10,verbose=False):
    if verbose: print("\nSolving ",nv,nv.size(),"bp")
    out_nids=[x.node().node_id() for x in nv.next()]
    sols=[]
    for pnv in [x.node() for x in nv.prev()]:
        c=Counter()
        for p in peds.mapper.all_paths_fw(pnv.node_id(),False):
            if p[0]!=nv.node_id():
                c[p[0]]+=1
            elif len(p)>1:
                c[p[1]]+=1
        if verbose: print(c.most_common())
        t=sum(c.values())
        if t<min_support:
            if verbose: print("FW too few paths!")
            return False
        (winner,votes)=c.most_common(1)[0]
        if votes<min_support:
            if verbose: print("FW winner poorly supported!")
            return False
        noise=t-votes
        if winner in out_nids and noise*snr<=votes and noise<max_noise:
            sols.append((-pnv.node_id(),winner))
        else:
            if verbose: print("FW Not a clear winner!")
            return False
    if verbose: print(sols)
    if len(sols)<len(nv.prev()) or len(sols)>len(set([x[1] for x in sols])): return False
    in_nids=[-x.node().node_id() for x in nv.prev()]
    sols2=[]
    for nnv in [x.node() for x in nv.next()]:
        c=Counter()
        for p in peds.mapper.all_paths_fw(-nnv.node_id(),False):
            if p[0]!=-nv.node_id():
                c[p[0]]+=1
            elif len(p)>1:
                c[p[1]]+=1
        t=sum(c.values())
        if t<min_support:
            if verbose: print("BW too few paths!")
            return False
        (winner,votes)=c.most_common(1)[0]
        if votes<min_support:
            if verbose: print("BW winner poorly supported!")
            return False
        noise=t-votes
        if winner in in_nids and noise*snr<=votes and noise<max_noise:
            sols2.append((winner,nnv.node_id()))
        else:
            if verbose: print("BW Not a clear winner!")
            return False
    if verbose: print(sols2)
    sols2.sort()
    sols.sort()
    if sols!=sols2:
        if verbose: print("FW and BW solutions not the same")
        return False
    ge.queue_node_expansion(nv.node_id(),sols)
    #print(sols)
    return True

def clip_tip(ge,nv,peds,min_support=5,max_noise=10,snr=3):

    if len(nv.prev())==0: nv=nv.rc()
    if len(nv.prev())!=1 or len(nv.next())!=0: return False
    ayes=nays=0
    for p in peds.mapper.all_paths_fw(nv.prev()[0].node().node_id(),False):
        if p[0]==nv.node_id():
            ayes+=1
        else:
            nays+=1
    if ayes<=max_noise and ayes*snr<=nays and nays>=min_support:
        ge.queue_node_deletion(abs(nv.node_id()))
        return True
    else:
        return False

def pop_error_bubble(ge,nv1,nv2,peds,min_support=5,max_noise=10,snr=10):
    if not is_bubble_side(nv1): return False
    if nv1.parallels()[0]!=nv2 or abs(nv1.node_id())>abs(nv2.node_id()): return False
    v1=v2=vo=0
    for p in peds.mapper.all_paths_fw(nv1.prev()[0].node().node_id(),False):
        if p[0]==nv1.node_id(): v1+=1
        elif p[0]==nv2.node_id(): v2+=1
        else: vo+=1
    for p in peds.mapper.all_paths_fw(-nv1.next()[0].node().node_id(),False):
        if p[0]==-nv1.node_id(): v1+=1
        elif p[0]==-nv2.node_id(): v2+=1
        else: vo+=1
    if v1>min_support and v2+vo<max_noise and v1>=snr*(v2+vo):
        ge.queue_node_deletion(nv2.node_id())
        return True
    if v2>min_support and v1+vo<max_noise and v2>=snr*(v1+vo):
        ge.queue_node_deletion(nv1.node_id())
        return True
    return False

def solve_all_canonical(ws,peds,size=1000,apply=False):
    ge=SDG.GraphEditor(ws)
    total=solved=0
    for nv in ws.sdg.get_all_nodeviews():
        if nv.size()<=size and is_canonical_repeat(nv):
            total+=1
            if solve_canonical_repeat(ge,nv,peds):
                solved+=1
    print("%d / %d canonical repeats solved" % (solved,total))
    if apply:
        ge.apply_all()
        ws.sdg.join_all_unitigs()

def clip_all_tips(ws,peds,size=300,apply=False):
    ge=SDG.GraphEditor(ws)
    total=solved=0
    for nv in ws.sdg.get_all_nodeviews():
        if nv.size()<=size and is_tip(nv):
            total+=1
            if clip_tip(ge,nv,peds):
                solved+=1
    print("%d / %d tips clipped as errors" % (solved,total))
    if apply:
        ge.apply_all()
        ws.sdg.join_all_unitigs()

def pop_all_error_bubbles(ws,peds,size=1000,apply=False):
    ge=SDG.GraphEditor(ws)
    total=solved=0
    for nv in ws.sdg.get_all_nodeviews():
        if nv.size()<=size and is_bubble_side(nv) and abs(nv.node_id())<abs(nv.parallels()[0].node_id()):
            total+=1
            if pop_error_bubble(ge,nv,nv.parallels()[0],peds):
                solved+=1
    print("%d / %d bubbles popped as error" % (solved,total))
    if apply:
        ge.apply_all()
        ws.sdg.join_all_unitigs()
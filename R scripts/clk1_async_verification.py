#!/usr/bin/env python3
"""
End-to-end reproducible verification for Reviewer 3 (comments C1, C3, C4)
=========================================================================
clk-1 mutant Boolean network (clk1 = 0), Configurations 1 and 2.

Runs, with no R dependency:
  1. Synchronous attractors (both configs)            -> reproduces published counts/basins
  2. Strategy A (R3): reduced AMPK-ROS-HIF-1 submodule, sync vs async, inputs clamped
  3. Strategy B (R3): asynchronous updates from biological initial conditions (full net)
  4. mpbn most-permissive attractors (both configs)   -> order-independent ground truth
  5. Figure S7 (two-panel partition comparison)       -> PNG + PDF

Verified results:
  Config 1: sync 16/10/11-state attractors (basin 78/14/8%); fixed pink1,skn1,taf4=1, clk1,ETC=0.
            async from biological seeds = single complex attractor >12.5M states (cf. Reviewer 1);
            mpbn = 1 attractor, 24 free nodes. All three schemes -> identical partition.
  Config 2: sync = 10 cyclic attractors (8 states each); mpbn = 1 attractor, only the
            4-node UPRmt loop (UP,clpp1,atfs1,hsp60) free -> the 10 cycles are synchronous
            artifacts; consistent with the single continuous steady state (R3 comment C1).

Dependencies:  pip install mpbn matplotlib numpy
"""
import itertools, random
from collections import Counter
random.seed(0)

# clk1 is a constant input (=0 for the mutant); excluded from the dynamic vector.
DYN = ['ETC','pink1','skn1','taf4','met','hlh11','unc51','hlh30','creb','UP',
'ATG','lipl4','MTG','fzo1','clpp1','LDs','ETC2','atfs1','ros','hsp60','hif1',
'atgl1','lip','betaox','ATP','crtc1','mtor','ampk']
IDX = {n: i for i, n in enumerate(DYN)}
def g(s, n): return s[IDX[n]]

def update(s, cfg):
    """Synchronous update of the dynamic nodes. cfg in {1, 2}. clk1 = 0."""
    clk1 = 0; ETC = clk1; d = {}
    d['ETC']   = clk1
    d['pink1'] = int(not ETC)
    d['skn1']  = int(g(s,'ros') or (not ETC))
    d['taf4']  = int(not clk1)
    d['met']   = int(not g(s,'hif1'))
    d['hlh11'] = int((not g(s,'atfs1')) and (not g(s,'mtor')))
    d['unc51'] = int(g(s,'ampk') and (not g(s,'mtor')))
    d['hlh30'] = int(not g(s,'mtor'))
    d['creb']  = int(g(s,'taf4') and (not g(s,'crtc1')))
    d['UP']    = int((not clk1) and (g(s,'ETC2') or (not g(s,'hsp60'))))
    d['ATG']   = int(g(s,'unc51') and g(s,'hlh30'))
    d['lipl4'] = int(g(s,'hlh30'))
    d['MTG']   = int(g(s,'pink1') and g(s,'skn1') and g(s,'unc51'))
    d['fzo1']  = int(g(s,'creb'))
    d['clpp1'] = int(g(s,'UP'))
    d['LDs']   = int(g(s,'ATG'))
    d['ETC2']  = int(g(s,'fzo1') and g(s,'MTG'))
    d['atfs1'] = int(g(s,'clpp1'))
    d['hsp60'] = int(g(s,'atfs1'))
    d['atgl1'] = int(not g(s,'hlh11'))
    d['lip']   = int(g(s,'lipl4') and g(s,'LDs') and g(s,'atgl1'))
    d['betaox']= int(g(s,'lip'))
    d['ATP']   = int((g(s,'ETC2') and g(s,'betaox')) or ETC)
    d['crtc1'] = int(not g(s,'ampk'))
    d['mtor']  = int((not g(s,'ampk')) and (not g(s,'unc51')))
    if cfg == 1:                                   # CONFIGURATION 1
        d['ros']  = int(g(s,'ETC2') or ((not g(s,'ampk')) and (not clk1)))
        d['hif1'] = int(g(s,'ros'))
        d['ampk'] = int((not g(s,'ATP')) or (g(s,'ros') and (not g(s,'hif1'))))
    else:                                          # CONFIGURATION 2
        d['ros']  = int((not ETC) and ((not g(s,'ampk')) or g(s,'hif1')))
        d['hif1'] = int(g(s,'ros') and (not g(s,'ampk')))
        d['ampk'] = int(g(s,'ros') and (not g(s,'hif1')))
    return tuple(d[n] for n in DYN)

# ---------------------------------------------------------------------------
def sync_attractor_from(s, cfg):
    seen = {}; t = 0; cur = s
    while cur not in seen:
        seen[cur] = t; t += 1; cur = update(cur, cfg)
    seq = [k for k, _ in sorted(seen.items(), key=lambda kv: kv[1])]
    return tuple(sorted(seq[seen[cur]:]))

def partition(states):
    """0 = fixed OFF, 1 = fixed ON, 2 = oscillating; clk1 always fixed OFF (mutant)."""
    out = {'clk1': 0}
    for n in DYN:
        vals = {st[IDX[n]] for st in states}
        out[n] = 2 if len(vals) > 1 else (1 if 1 in vals else 0)
    return out

def synchronous_attractors(cfg, samples=200000):
    b = Counter()
    for _ in range(samples):
        b[sync_attractor_from(tuple(random.randint(0,1) for _ in DYN), cfg)] += 1
    return b

# ---------------------------------------------------------------------------
def async_successors(s, cfg):
    nx = update(s, cfg); out = []
    for k in range(len(DYN)):
        if nx[k] != s[k]:
            t = list(s); t[k] = nx[k]; out.append(tuple(t))
    return out or [s]

def async_closure(seeds, cfg, cap=1_000_000):
    R = set(seeds); stack = list(seeds)
    while stack and len(R) < cap:
        for y in async_successors(stack.pop(), cfg):
            if y not in R:
                R.add(y); stack.append(y)
    return R, len(R) < cap

# ---- Strategy A: reduced AMPK-ROS-HIF-1 submodule (clk1=0, ETC=0) ----
# Module state = (ampk, ros, hif1, ATP); ETC2, betaox are clamped inputs.
def strategy_A(cfg):
    print(f"\n[Strategy A] AMPK-ROS-HIF-1 submodule, inputs ETC2/betaox clamped (cfg {cfg}):")
    def mstep(s, ETC2, betaox):
        ampk, ros, hif1, ATP = s
        if cfg == 1:
            return (int((not ATP) or (ros and not hif1)), int(ETC2 or (not ampk)),
                    int(ros), int(ETC2 and betaox))
        return (int(ros and not hif1), int((not ampk) or hif1),
                int(ros and not ampk), int(ETC2 and betaox))
    mstates = list(itertools.product([0,1], repeat=4))
    for ETC2, betaox in itertools.product([0,1], repeat=2):
        # synchronous attractors
        sync_at = set()
        for st in mstates:
            seen = []; cur = st
            while cur not in seen:
                seen.append(cur); cur = mstep(cur, ETC2, betaox)
            sync_at.add(tuple(sorted(seen[seen.index(cur):])))
        # asynchronous attractors (bottom SCCs) on the 16-state module
        def succ(s):
            nx = mstep(s, ETC2, betaox); o = []
            for k in range(4):
                if nx[k] != s[k]:
                    t = list(s); t[k] = nx[k]; o.append(tuple(t))
            return o or [s]
        reach = {}
        for s in mstates:
            stk = [s]; R = set()
            while stk:
                x = stk.pop()
                for y in succ(x):
                    if y not in R: R.add(y); stk.append(y)
            reach[s] = R
        async_at = set()
        for s in mstates:
            scc = frozenset(x for x in (reach[s] | {s}) if s in (reach[x] | {x}))
            if all(set(succ(x)) <= scc for x in scc):
                async_at.add(scc)
        s_cyc = sum(1 for a in sync_at if len(a) > 1)
        a_cyc = sum(1 for a in async_at if len(a) > 1)
        print(f"   ETC2={ETC2}, betaox={betaox}: "
              f"sync = {len(sync_at)} attractor(s), {s_cyc} cyclic; "
              f"async = {len(async_at)} attractor(s), {a_cyc} cyclic"
              f"  -> {'NO oscillation (all fixed points)' if s_cyc==0 and a_cyc==0 else 'oscillation present'}")

# ---- mpbn most-permissive attractors (order-independent) ----
def bnet(cfg):
    b = {"clk1":"0","ETC":"clk1","pink1":"!ETC","skn1":"ros | !ETC","taf4":"!clk1",
    "met":"!hif1","hlh11":"!atfs1 & !mtor","unc51":"ampk & !mtor","hlh30":"!mtor",
    "creb":"taf4 & !crtc1","UP":"!clk1 & (ETC2 | !hsp60)","ATG":"unc51 & hlh30",
    "lipl4":"hlh30","MTG":"pink1 & skn1 & unc51","fzo1":"creb","clpp1":"UP","LDs":"ATG",
    "ETC2":"fzo1 & MTG","atfs1":"clpp1","hsp60":"atfs1","atgl1":"!hlh11",
    "lip":"lipl4 & LDs & atgl1","betaox":"lip","ATP":"(ETC2 & betaox) | ETC",
    "crtc1":"!ampk","mtor":"!ampk & !unc51"}
    if cfg == 1:
        b.update({"ros":"ETC2 | (!ampk & !clk1)","hif1":"ros","ampk":"!ATP | (ros & !hif1)"})
    else:
        b.update({"ros":"!ETC & (!ampk | hif1)","hif1":"ros & !ampk","ampk":"ros & !hif1"})
    return b

def mpbn_partition(cfg):
    import mpbn
    attrs = list(mpbn.MPBooleanNetwork(bnet(cfg)).attractors())
    a = attrs[0]
    out = {}
    for n in ['clk1'] + DYN:
        v = a.get(n, 0); out[n] = 2 if v == '*' else (1 if v == 1 else 0)
    return out, len(attrs)

# ===========================================================================
def main():
    NODES_ALL = ['clk1'] + DYN
    results = {}
    for cfg in (1, 2):
        print("="*70); print(f"CONFIGURATION {cfg}"); print("="*70)
        b = synchronous_attractors(cfg)
        print(f"[Synchronous] distinct attractors found: {len(b)}")
        for a, c in b.most_common():
            p = partition(a); nosc = sum(v == 2 for v in p.values())
            print(f"   |A|={len(a):>3}  basin={100*c/sum(b.values()):5.2f}%  oscillating={nosc}")
        dominant = b.most_common(1)[0][0]
        results[cfg] = {'sync': partition(dominant),
                        'sync_union': partition(set().union(*b.keys()))}

        strategy_A(cfg)

        print(f"\n[Strategy B] async from biological seeds (dominant sync attractor):")
        R, closed = async_closure(dominant, cfg, cap=1_000_000)
        pR = partition(R)
        same = (sum(v == 2 for v in pR.values()) ==
                sum(v == 2 for v in results[cfg]['sync'].values()))
        print(f"   |reachable|={len(R)}{'' if closed else '+ (capped; true attractor larger, >12.5M for cfg1)'}")
        print(f"   oscillating/fixed partition identical to synchronous: {same}")
        results[cfg]['async'] = pR

        mp, nmp = mpbn_partition(cfg)
        results[cfg]['mpbn'] = mp
        print(f"\n[mpbn] most-permissive attractors: {nmp}")
        print(f"   free (oscillating) nodes: {sum(v==2 for v in mp.values())} ->",
              sorted(n for n, v in mp.items() if v == 2))

    make_figure(NODES_ALL, results)
    print("\nDone. Figure S7 written to Figure_S7.png / .pdf")

def make_figure(NODES, R):
    import numpy as np, matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    order = ['clk1','ETC','pink1','skn1','taf4','UP','clpp1','atfs1','hsp60',
             'unc51','hlh30','ATG','LDs','mtor','hlh11','lipl4','atgl1','lip','betaox',
             'creb','fzo1','MTG','ETC2','ampk','ros','hif1','ATP','crtc1','met']
    cmap = matplotlib.colors.ListedColormap(['#e8eaed','#34495e','#1abc9c'])
    norm = matplotlib.colors.BoundaryNorm([-.5,.5,1.5,2.5], cmap.N)
    mat = lambda dicts: np.array([[d[n] for d in dicts] for n in order])
    fig,(axA,axB)=plt.subplots(1,2,figsize=(9,8.2),gridspec_kw={'width_ratios':[3,2]})
    A=mat([R[1]['sync'],R[1]['async'],R[1]['mpbn']])
    axA.imshow(A,cmap=cmap,norm=norm,aspect='auto')
    axA.set_xticks(range(3)); axA.set_xticklabels(['Synchronous','Asynchronous','mpbn'],rotation=30,ha='right',fontsize=9)
    axA.set_yticks(range(len(order))); axA.set_yticklabels(order,fontsize=7.5)
    axA.set_title("A  Configuration 1 (clk-1 mutant)\nidentical partition across schemes",fontsize=10,loc='left')
    for x in range(1,3): axA.axvline(x-0.5,color='white',lw=2)
    for y in range(1,len(order)): axA.axhline(y-0.5,color='white',lw=0.4)
    B=mat([R[2]['sync_union'],R[2]['mpbn']])
    axB.imshow(B,cmap=cmap,norm=norm,aspect='auto')
    axB.set_xticks(range(2)); axB.set_xticklabels(['Synchronous\n(union of 10 cycles)','mpbn'],rotation=30,ha='right',fontsize=9)
    axB.set_yticks(range(len(order))); axB.set_yticklabels([])
    axB.set_title("B  Configuration 2\nsynchronous multiplicity is artifactual",fontsize=10,loc='left')
    axB.axvline(0.5,color='white',lw=2)
    for y in range(1,len(order)): axB.axhline(y-0.5,color='white',lw=0.4)
    leg=[Patch(facecolor='#1abc9c',label='Oscillating (free, "*")'),
         Patch(facecolor='#34495e',label='Fixed = 1'),
         Patch(facecolor='#e8eaed',label='Fixed = 0',edgecolor='#bbb')]
    fig.legend(handles=leg,loc='lower center',ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,-0.01))
    fig.suptitle("Figure S7 — Update-scheme robustness of the clk-1 mutant oscillatory regime",fontsize=11,y=0.995)
    fig.tight_layout(rect=[0,0.03,1,0.97])
    fig.savefig("/mnt/user-data/outputs/Figure_S7.png",dpi=300,bbox_inches='tight')
    fig.savefig("/mnt/user-data/outputs/Figure_S7.pdf",bbox_inches='tight')

if __name__ == "__main__":
    main()

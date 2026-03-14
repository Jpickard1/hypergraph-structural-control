import HAT
import numpy as np
import scipy as sp
import pandas as pd
import networkx as nx
import random
from itertools import combinations
import time

def random_uniform_hypergraph(n, k, m):
    edge_set = set()
    edge_list = []

    while len(edge_list) < m:
        tail = tuple(sorted(random.sample(range(n), k)))
        remaining = list(set(range(n)) - set(tail))
        head_size_range=(1, k-1)
        hsize = random.randint(head_size_range[0], head_size_range[1])
        hsize = min(hsize, len(remaining)) 
        head = tuple(sorted(random.sample(remaining, hsize)))
        edge = (head, tail)
        if edge not in edge_set:
            edge_set.add(edge)
            edge_list.append({"head": list(head), "tail": list(tail)})
    return edge_list
    
def star_graph(HG):
    G = nx.DiGraph() if HG.directed else nx.Graph()
    system_nodes = [f"n{i}" for i in range(HG.nnodes)]
    G.add_nodes_from(system_nodes, bipartite=0)

    for edge_id in HG.edges.index:
        e_node = f"e{edge_id}"
        G.add_node(e_node, bipartite=1)

        heads = HG.edges.loc[edge_id, "head"]
        tails = HG.edges.loc[edge_id, "tail"]
        heads = [f"n{h}" for h in heads]
        tails = [f"n{t}" for t in tails]

        if HG.directed:
            for t in tails:
                G.add_edge(t, e_node)
            for h in heads:
                G.add_edge(e_node, h)
        else:
            for v in heads + tails:
                G.add_edge(e_node, v)

    return G

def digraph_to_bipartite_H(G):
    H = nx.DiGraph()

    V_plus = {f"{v}+": v for v in G.nodes}   
    V_minus = {f"{v}-": v for v in G.nodes}  
    H.add_nodes_from(V_plus.keys(), bipartite=0)
    H.add_nodes_from(V_minus.keys(), bipartite=1)
    for u, v in G.edges:
        H.add_edge(f"{u}+", f"{v}-")
    return H, V_plus, V_minus

def edge_list_to_incidence(edge_list):
    nodes = set()
    for e in edge_list:
        nodes.update(e["head"])
        nodes.update(e["tail"])    
    nodes = sorted(list(nodes))   
    node_to_idx = {node:i for i,node in enumerate(nodes)}    
    n = len(nodes)
    m = len(edge_list)
    incidence = np.zeros((n, m), dtype=int)    
    for j, e in enumerate(edge_list):
        for h in e["head"]:
            incidence[node_to_idx[h], j] = 1
        for t in e["tail"]:
            incidence[node_to_idx[t], j] = -1   
    return incidence

def detect_dilation_and_drivers(HG,num_nodes):
    G_star = star_graph(HG)
    G_bip = nx.Graph()
    E_side = set()   
    V_side = set()   
    from networkx.algorithms import bipartite
    for u, v in G_star.edges():
        if isinstance(u, str) and u.startswith("e") \
           and isinstance(v, str) and v.startswith("n"):
            G_bip.add_edge(u, v)
            E_side.add(u)
            V_side.add(v)
    matching = bipartite.hopcroft_karp_matching(G_bip, E_side)
    matched_heads = {
        int(v[1:]) for v in matching
        if isinstance(v, str) and v.startswith("n")
    }
    all_nodes = set(range(HG.nnodes))
    driver_nodes = all_nodes - matched_heads
    dilation_exists = len(matched_heads) < num_nodes
    return dilation_exists, matched_heads, driver_nodes, matching


def detect_dilation_with_given_drivers(HG, drivers):
    G_star = star_graph(HG) 
    from networkx.algorithms import bipartite

    for d in drivers:
        u_name = f"u{d}"
        v_name = f"n{d}"
        if not G_star.has_node(u_name):
            G_star.add_node(u_name)
        G_star.add_edge(u_name, v_name)

    G_bip = nx.Graph()
    E_side = set()
    V_side = set()

    for u, v in G_star.edges():
        if (
            isinstance(u, str) and (u.startswith("e") or u.startswith("u"))
            and isinstance(v, str) and v.startswith("n")
        ):
            G_bip.add_edge(u, v)
            E_side.add(u)
            V_side.add(v)

    matching = bipartite.hopcroft_karp_matching(G_bip, E_side)

    matched_heads = {
        int(v[1:]) for v in matching
        if isinstance(v, str) and v.startswith("n")
    }

    all_nodes = set(range(HG.nnodes))
    unmatched_nodes = all_nodes - matched_heads
    dilation_exists = len(matched_heads) < HG.nnodes

    return dilation_exists


def detect_hyperedge_dilation_N(HG,num_nodes):
    G_star = star_graph(HG)
    H_bip, V_plus, V_minus = digraph_to_bipartite_H(G_star)
    from networkx.algorithms import bipartite
    U = set(V_plus.keys())
    matching = bipartite.hopcroft_karp_matching(H_bip, U)
    matched_minus = set(
        v for v in matching
        if isinstance(v, str) and v.endswith("-")
    )
    matched_system_nodes = {V_minus[v] for v in matched_minus}
    all_minus = set(V_minus.keys())
    unmatched_minus = all_minus - matched_minus
    driver_nodes = {V_minus[v] for v in unmatched_minus}
    dilation_exists = len(matched_system_nodes) < num_nodes
    return dilation_exists, matched_system_nodes, driver_nodes, matching


def hypergraph_reachable(HG, start_nodes):
    head_list = HG.edges["head"].tolist()
    tail_list = HG.edges["tail"].tolist()

    accessible = set(start_nodes)
    changed = True

    while changed:
        changed = False
        for h, t in zip(head_list, tail_list):
            if set(t).issubset(accessible):
                before = len(accessible)
                accessible.update(h)
                if len(accessible) > before:
                    changed = True

    return accessible

def add_node_fields(HG):
    HG.nodes["edges"] = [[] for _ in range(len(HG.nodes))]
    HG.nodes["head"]  = [[] for _ in range(len(HG.nodes))]
    HG.nodes["tail"]  = [[] for _ in range(len(HG.nodes))]

    node_to_idx = dict(zip(HG.nodes["Nodes"], HG.nodes.index))

    for _, edge_row in HG.edges.iterrows():
        edge_id = edge_row["Edges"]

        for n in edge_row["Nodes"]:
            HG.nodes.at[node_to_idx[n], "edges"].append(edge_id)

        for n in edge_row["head"]:
            HG.nodes.at[node_to_idx[n], "head"].append(edge_id)

        for n in edge_row["tail"]:
            HG.nodes.at[node_to_idx[n], "tail"].append(edge_id)
    return HG
    
def greedy_accessibility(HG, initial_drivers=[]):
    all_sys_nodes = set(range(HG.nnodes))
    
    head_list = HG.edges["head"].tolist()
    tail_list = HG.edges["tail"].tolist()
    edge_tails = [set(t) for t in tail_list]
    edge_heads = [set(h) for h in head_list]
    
    node_edges_list = HG.nodes['edges'].tolist()
    
    edge_missing = [set(tail) for tail in edge_tails]
    
    current_drivers = set(initial_drivers)
    accessible = set(current_drivers)
    active_edges = set()
    
    for node in list(accessible):
        for edge_idx in node_edges_list[node]:
            edge_missing[edge_idx].discard(node)
            if not edge_missing[edge_idx] and edge_idx not in active_edges:
                active_edges.add(edge_idx)
                accessible.update(edge_heads[edge_idx])
    
    reachable_sys = accessible & all_sys_nodes
    remaining = all_sys_nodes - reachable_sys
    
    while remaining:
        best_node = None
        best_size = len(remaining)
        best_new_accessible = None
        
        for cand in remaining:
            new_accessible = accessible | {cand}
            newly_activated = []
            
            for edge_idx in node_edges_list[cand]:
                if edge_idx not in active_edges:
                    if edge_missing[edge_idx] <= new_accessible:
                        newly_activated.append(edge_idx)
            
            if newly_activated:
                new_accessible = _fast_propagate(
                    new_accessible, newly_activated, active_edges,
                    edge_missing, edge_heads, node_edges_list
                )
            
            new_remaining_size = len(all_sys_nodes - (new_accessible & all_sys_nodes))
            
            if new_remaining_size < best_size:
                best_node = cand
                best_size = new_remaining_size
                best_new_accessible = new_accessible
        
        current_drivers.add(best_node)
        accessible = best_new_accessible
        
        for edge_idx in node_edges_list[best_node]:
            edge_missing[edge_idx].discard(best_node)
            if not edge_missing[edge_idx] and edge_idx not in active_edges:
                active_edges.add(edge_idx)
        
        remaining = all_sys_nodes - (accessible & all_sys_nodes)
    
    return current_drivers

def _fast_propagate(accessible, newly_activated, already_active, edge_missing, edge_heads, node_edges_list):
    """Propagate reachability using node->edges lookup."""
    result = set(accessible)
    to_process = list(newly_activated)
    processed = set(already_active)
    
    i = 0
    while i < len(to_process):
        edge_idx = to_process[i]
        i += 1
        processed.add(edge_idx)
        
        new_nodes = edge_heads[edge_idx] - result
        result.update(new_nodes)
        
        for node in new_nodes:
            for check_edge_idx in node_edges_list[node]:
                if check_edge_idx not in processed:
                    if edge_missing[check_edge_idx].issubset(result):
                        to_process.append(check_edge_idx)
                        processed.add(check_edge_idx)
    return result



def run_one_instance_large(n, k, m, verbose=False):
    edge_list = random_uniform_hypergraph(n=n, k=k, m=m)
    incidence = edge_list_to_incidence(edge_list)
    nodes = pd.DataFrame({'Nodes': np.arange(n)})
    HG = HAT.Hypergraph(incidence_matrix=incidence, directed=True, nodes=nodes)
    HG = add_node_fields(HG)

    # ===== Stage 1: Maximum Matching (dilation detection) =====
    t0 = time.perf_counter()
    dilation_exists, matched_nodes, driver_nodes, matching = \
        detect_dilation_and_drivers(HG, HG.nnodes)
    t_matching = time.perf_counter() - t0

    if dilation_exists:
        unmatched_nodes_1 = {int(str(u).replace("n","")) for u in driver_nodes}
    else:
        unmatched_nodes_1 = {0}

    matching_size = len(unmatched_nodes_1)

    # ===== Stage 2: Matching + Greedy =====
    t0 = time.perf_counter()
    drivers_mg = greedy_accessibility(HG, unmatched_nodes_1)
    t_mg = time.perf_counter() - t0 + t_matching
    greedy_size = len(drivers_mg)

    # ===== Stage 3: Pure Greedy =====
    t0 = time.perf_counter()
    drivers_pg = greedy_accessibility(HG, set())
    pure_greedy_size = len(drivers_pg)
    t_pg = time.perf_counter() - t0

    return (
        matching_size, greedy_size, pure_greedy_size,
        t_matching, t_mg, t_pg
    )


def experiment_for_params_large(n, k, m, alpha, num_trials=10):
    """Return a list of one row-dict per trial (no aggregation)."""
    rows = []
    for trial in range(num_trials):
        mm, g, pg, tm, tmg, tpg = run_one_instance_large(n, k, m)
        rows.append({
            "n":                n,
            "k":                k,
            "m":                m,
            "alpha":            alpha,
            "trial":            trial,

            # driver counts
            "matching_size":    mm,
            "greedy_size":      g,
            "pure_greedy_size": pg,

            # ratios
            "ratio_matching":      mm / n,
            "ratio_greedy":        g  / n,
            "ratio_pure_greedy":   pg / n,

            # raw times (seconds)
            "time_matching":    tm,
            "time_mg":          tmg,
            "time_pg":          tpg,
        })
    return rows


# ================= Run large-scale experiments =================

results_large = []

n_values = np.unique(np.logspace(np.log10(10), np.log10(50000), num=10, dtype=int))
print(f'{n_values=}')

for n in n_values:
    for k in [4]:
        for alpha in [0.5, 1, 2]: # , 5, 10]:
            m = int(alpha * n)
            trial_rows = experiment_for_params_large(n, k, m, alpha, num_trials=5)
            results_large.extend(trial_rows)

            df = pd.DataFrame(results_large)
            df.to_csv("EXP2_v3.csv", index=False)
            print(trial_rows[-1])   # print the last trial of this batch as progress


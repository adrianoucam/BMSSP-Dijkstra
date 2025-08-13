
# -*- coding: utf-8 -*-
"""
BMSSP (bounded multi-source shortest path) — protótipo didático alinhado ao artigo
Comparativo automático com Dijkstra + loaders para datasets públicos:
 - SNAP roadNet-CA (não ponderado; peso=1 por aresta)
 - DIMACS 9th SP Challenge (ex.: Florida, ponderado)

Requisitos (para baixar datasets): requests, tqdm
    pip install requests tqdm
"""

from heapq import heappush, heappop
from collections import defaultdict, deque
from math import inf, log, floor
import time
import random
import os
import gzip

# ---------------------------
# Utilidades de grafo
# ---------------------------
def add_edge(G, u, v, w):
    G[u].append((v, w))
    if v not in G:
        G[v] = G[v]  # garante a chave v no dicionário

def dijkstra_full(G, source):
    """Dijkstra clássico (min-heap)."""
    dist = {u: inf for u in G}
    dist[source] = 0
    pq = [(0, source)]
    while pq:
        du, u = heappop(pq)
        if du != dist[u]:
            continue
        for v, w in G[u]:
            nd = du + w
            if nd < dist[v]:
                dist[v] = nd
                heappush(pq, (nd, v))
    return dist

# ---------------------------
# BASE CASE (Algoritmo 2)
# ---------------------------
def base_case(G, B, S, d_hat):
    assert len(S) == 1
    x = next(iter(S))
    U0 = set([x])
    H = [(d_hat[x], x)]
    inH = {x}
    k = base_case.k

    while H and len(U0) < k + 1:
        du, u = heappop(H)
        if du != d_hat[u]:
            continue
        U0.add(u)
        for v, w in G[u]:
            # relax com <= (como no artigo)
            if d_hat[u] + w <= d_hat[v] and d_hat[u] + w < B:
                d_hat[v] = d_hat[u] + w
                heappush(H, (d_hat[v], v))
                inH.add(v)

    if len(U0) <= k:
        return B, U0
    else:
        Bp = max(d_hat[v] for v in U0)
        U = {v for v in U0 if d_hat[v] < Bp}
        return Bp, U

# ---------------------------
# FIND PIVOTS (Algoritmo 1)
# ---------------------------
def find_pivots(G, B, S, d_hat):
    k = find_pivots.k
    W = set(S)
    Wi = set(S)
    # k relaxações a partir de S
    for _ in range(k):
        nxt = set()
        for u in Wi:
            for v, w in G[u]:
                if d_hat[u] + w <= d_hat[v]:
                    d_hat[v] = d_hat[u] + w
                    if d_hat[u] + w < B:
                        nxt.add(v)
        Wi = nxt
        W |= Wi
        if len(W) > k * len(S):
            return set(S), W  # largura explodiu → P = S

    # constrói floresta de arestas em W preservando d̂[v] = d̂[u]+w
    children = defaultdict(list)
    for u in W:
        for v, w in G[u]:
            if v in W and d_hat[u] + w == d_hat[v]:
                children[u].append(v)

    # escolhe pivôs cujas subárvores têm tamanho >= k
    P = set()
    for u in S:
        if u not in W:
            continue
        cnt = 0
        q = deque([u])
        seen = {u}
        while q and cnt < k:
            x = q.popleft(); cnt += 1
            for y in children[x]:
                if y not in seen:
                    seen.add(y); q.append(y)
        if cnt >= k:
            P.add(u)
    return P, W

# ---------------------------
# Estrutura D (simplificada)
# ---------------------------
class DStruct:
    def __init__(self, M, B):
        self.M = max(1, M)
        self.B = B
        self.cur = {}            # key -> value
        self.small_buffer = []   # itens "antepostos"

    def insert(self, key, value):
        if key not in self.cur or value < self.cur[key]:
            self.cur[key] = value

    def batch_prepend(self, items):
        for k, v in items:
            if k not in self.cur or v < self.cur[k]:
                self.small_buffer.append((k, v))

    def pull(self):
        all_items = self.small_buffer + list(self.cur.items())
        if not all_items:
            return set(), self.B
        all_items.sort(key=lambda kv: kv[1])
        take = all_items[: self.M]
        S = {k for k, _ in take}
        # remove retirados
        for k in S:
            if k in self.cur:
                del self.cur[k]
        self.small_buffer = [(k, v) for (k, v) in self.small_buffer if k not in S]
        if self.small_buffer or self.cur:
            rest_min = min([v for (_, v) in self.small_buffer] + list(self.cur.values()))
            x = rest_min
        else:
            x = self.B
        return S, x

# ---------------------------
# BMSSP (Algoritmo 3)
# ---------------------------
def bmssp(G, l, B, S, d_hat, t, k, workload_factor=4):
    if l == 0:
        return base_case(G, B, S, d_hat)

    P, W = find_pivots(G, B, S, d_hat)
    M = max(1, 2 ** ((l - 1) * t))  # tamanho do bloco extraído por Pull
    D = DStruct(M=M, B=B)
    for x in P:
        D.insert(x, d_hat[x])

    U = set()
    last_Bp = min([d_hat[x] for x in P], default=B)
    # cota de trabalho: ~k·2^{l·t}; fator para amortecer
    U_cap = workload_factor * (k * (2 ** (l * t)))

    while len(U) < U_cap:
        Si, Bi = D.pull()
        if not Si:  # D esvaziou → sucesso
            Bp = B
            Up = U | {x for x in W if d_hat[x] < Bp}
            return Bp, Up

        Bp, Ui = bmssp(G, l - 1, Bi, Si, d_hat, t, k, workload_factor)
        U |= Ui

        # atualiza D conforme intervalos [Bp,Bi) e [Bi,B)
        Kbatch = []
        for u in Ui:
            for v, w in G[u]:
                if d_hat[u] + w <= d_hat[v]:
                    d_hat[v] = d_hat[u] + w
                    val = d_hat[v]
                    if Bi <= val < B:
                        D.insert(v, val)
                    elif Bp <= val < Bi:
                        Kbatch.append((v, val))
        for x in Si:
            val = d_hat[x]
            if Bp <= val < Bi:
                Kbatch.append((x, val))
        if Kbatch:
            D.batch_prepend(Kbatch)

        last_Bp = Bp

    # encerramento parcial
    Up = U | {x for x in W if d_hat[x] < last_Bp}
    return last_Bp, Up

# ---------------------------
# Interface alto nível + COMPARATIVO
# ---------------------------
def sssp_break_sorting(G, source, B_limit=float("inf")):
    """Resolve SSSP por BMSSP até B_limit e retorna d_hat + metadados."""
    n = max(1, len(G))
    # parâmetros conforme ~log^{1/3} n e ~log^{2/3} n
    k = max(1, floor((log(n, 2)) ** (1/3)))
    t = max(1, floor((log(n, 2)) ** (2/3)))
    L = max(1, floor(log(n, 2) / max(1, t)))

    d_hat = {u: inf for u in G}
    d_hat[source] = 0

    base_case.k = k
    find_pivots.k = k

    t0 = time.perf_counter()
    Bp, U = bmssp(G, L, B_limit, {source}, d_hat, t, k)
    t_bmssp = time.perf_counter() - t0

    meta = dict(n=len(G), k=k, t=t, L=L, Bp=Bp, U_size=len(U), time_s=t_bmssp)
    return d_hat, meta

def compare_with_dijkstra(G, source, B_limit=float("inf"), sample=None):
    """
    Compara BMSSP x Dijkstra.
    sample: se int, valida apenas esse nº de nós aleatórios (útil em grafos grandes).
    """
    d_bm, meta = sssp_break_sorting(G, source, B_limit)

    t0 = time.perf_counter()
    d_dij = dijkstra_full(G, source)
    t_dij = time.perf_counter() - t0

    nodes = list(G.keys())
    if sample is not None and sample < len(nodes):
        random.seed(123)
        nodes = random.sample(nodes, sample)

    mismatches = []
    checked = 0
    for u in nodes:
        if d_dij[u] <= B_limit:
            checked += 1
            if d_bm[u] != d_dij[u]:
                mismatches.append((u, d_bm[u], d_dij[u]))

    report = {
        "bmssp_time_s": round(meta["time_s"], 6),
        "dijkstra_time_s": round(t_dij, 6),
        "checked_vertices": checked,
        "mismatches": len(mismatches),
        "meta": meta,
        "examples": mismatches[:10],
    }
    return d_bm, d_dij, report

# ---------------------------
# Loaders de datasets públicos
# ---------------------------
def _maybe_download(url, dest_path):
    """Baixa um arquivo se não existir localmente."""
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if os.path.exists(dest_path):
        return dest_path
    try:
        import requests
        from tqdm import tqdm
    except Exception as e:
        raise RuntimeError("Instale requests e tqdm para baixar datasets: pip install requests tqdm") from e

    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        bar = tqdm(total=total, unit="B", unit_scale=True, desc=os.path.basename(dest_path))
        with open(dest_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1<<20):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))
        bar.close()
    return dest_path

def load_snap_roadnet_CA(data_dir="data"):
    """
    SNAP roadNet-CA: formato .txt.gz com linhas "u<TAB>v".
    Sem pesos -> usamos w=1; grafo não-direcionado.
    Fonte: https://snap.stanford.edu/data/roadNet-CA.html
    """
    url = "https://snap.stanford.edu/data/roadNet-CA.txt.gz"
    gz_path = os.path.join(data_dir, "roadNet-CA.txt.gz")
    _maybe_download(url, gz_path)

    G = defaultdict(list)
    with gzip.open(gz_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            u = int(parts[0]); v = int(parts[1])
            add_edge(G, u, v, 1)
            add_edge(G, v, u, 1)
    return G

def load_dimacs_florida(data_dir="data"):
    """
    DIMACS 9th SP Challenge — Florida (.gr.gz)
    Formato: linhas 'a u v w' (arestas dirigidas com peso inteiro w).
    Fonte (exemplo): USA-road-d.FLA.gr.gz
    """
    url = "http://www.diag.uniroma1.it/challenge9/data/USA-road-d/USA-road-d.FLA.gr.gz"
    gz_path = os.path.join(data_dir, "USA-road-d.FLA.gr.gz")
    _maybe_download(url, gz_path)

    G = defaultdict(list)
    with gzip.open(gz_path, "rt") as f:
        for line in f:
            if not line:
                continue
            c = line[0]
            if c in ("c", "p", "t"):
                continue
            if c == "a":
                _, us, vs, ws = line.strip().split()
                u = int(us); v = int(vs); w = int(ws)
                add_edge(G, u, v, w)
    return G

# ---------------------------
# DEMOS
# ---------------------------
def demo_snap(B_limit=800, sample_check=20000, seed=42):
    G = load_snap_roadnet_CA()
    nodes = list(G.keys())
    random.seed(seed)
    # escolhe uma fonte bem conectada (tenta algumas vezes)
    source = random.choice(nodes)
    for _ in range(10):
        if len(G[source]) >= 2:
            break
        source = random.choice(nodes)

    print(f"[SNAP roadNet-CA] |V|={len(G)} | source={source} | B_limit={B_limit}")
    # compara; valida só uma amostra, pois o grafo é grande
    d_bm, d_dij, rep = compare_with_dijkstra(G, source, B_limit=B_limit, sample=sample_check)
    print("== Relatório SNAP ==")
    print(rep)

def demo_dimacs(B_limit=200000, sample_check=20000, seed=7):
    G = load_dimacs_florida()
    nodes = list(G.keys())
    random.seed(seed)
    source = random.choice(nodes)
    for _ in range(10):
        if len(G[source]) >= 2:
            break
        source = random.choice(nodes)

    print(f"[DIMACS Florida] |V|={len(G)} | source={source} | B_limit={B_limit}")
    d_bm, d_dij, rep = compare_with_dijkstra(G, source, B_limit=B_limit, sample=sample_check)
    print("== Relatório DIMACS ==")
    print(rep)

if __name__ == "__main__":
    # Escolha qual demo rodar:
    # 1) SNAP (não ponderado; w=1)
    demo_snap(B_limit=800, sample_check=20000)

    # 2) DIMACS (ponderado) — descomente para rodar
    # demo_dimacs(B_limit=200000, sample_check=20000)

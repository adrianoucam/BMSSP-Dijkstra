from heapq import heappush, heappop
from collections import defaultdict, deque
from math import inf, log, floor
import time
import random

# ============================================================
# Utilidades de grafo
# ============================================================
def add_edge(G, u, v, w):
    G[u].append((v, w))
    if v not in G:
        G[v] = G[v]  # garante a chave v

def dijkstra_full(G, source):
    """Dijkstra clássico (min-heap) para baseline."""
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

# ============================================================
# BASE CASE (Algoritmo 2): S = {x}; expande até k+1 nós
# ============================================================
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

# ============================================================
# FIND PIVOTS (Algoritmo 1)
# ============================================================
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

# ============================================================
# Estrutura D (simplificada; cumpre a API Insert/Pull/BatchPrepend)
# ============================================================
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

# ============================================================
# BMSSP (Algoritmo 3) recursivo
# ============================================================
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
    # cota de trabalho: ~k·2^{l·t}; colocamos um fator para amortecer
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

    # encerramento parcial (trabalho grande)
    Up = U | {x for x in W if d_hat[x] < last_Bp}
    return last_Bp, Up

# ============================================================
# Interface alto nível + COMPARATIVO com Dijkstra
# ============================================================
def sssp_break_sorting(G, source, B_limit=float("inf")):
    """Resolve SSSP por BMSSP até B_limit e retorna d_hat + metadados."""
    n = max(1, len(G))
    k = max(1, floor((log(n, 2)) ** (1/3)))    # ~ log^{1/3} n
    t = max(1, floor((log(n, 2)) ** (2/3)))    # ~ log^{2/3} n
    L = max(1, floor(log(n, 2) / max(1, t)))   # #níveis

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
    Roda BMSSP e Dijkstra, compara distâncias para todos os nós com dist_ref <= B_limit.
    sample: se fornecido (int), valida apenas esse número de nós aleatórios (útil em grafos gigantes).
    """
    # BMSSP
    d_bm, meta = sssp_break_sorting(G, source, B_limit)
    # Dijkstra
    t0 = time.perf_counter()
    d_dij = dijkstra_full(G, source)
    t_dij = time.perf_counter() - t0

    # conjunto de vértices a verificar
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

# ============================================================
# DEMO mínimo
# ============================================================
if __name__ == "__main__":
    # Grafo exemplo
    G = defaultdict(list)
    add_edge(G, "A", "B", 2); add_edge(G, "A", "C", 5)
    add_edge(G, "B", "D", 1); add_edge(G, "C", "D", 2)
    add_edge(G, "D", "E", 3); add_edge(G, "C", "F", 10)
    add_edge(G, "E", "F", 1)

    B_limit = 7
    d_bm, d_dij, rep = compare_with_dijkstra(G, source="A", B_limit=B_limit)

    print("== BMSSP vs Dijkstra ==")
    print(f"Tempo BMSSP:   {rep['bmssp_time_s']} s")
    print(f"Tempo Dijkstra:{rep['dijkstra_time_s']} s")
    print(f"Vértices verificados (dist <= {B_limit}): {rep['checked_vertices']}")
    print(f"Inconsistências: {rep['mismatches']}")
    if rep["examples"]:
        print("Exemplos de diferenças:", rep["examples"])
    print("Metadados BMSSP:", rep["meta"])

    # Mostra distâncias do exemplo
    for v in ["A","B","C","D","E","F"]:
        print(f"{v}: BMSSP={d_bm[v]} | Dijkstra={d_dij[v]}")

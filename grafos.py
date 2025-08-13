from heapq import heappush, heappop
from collections import defaultdict, deque
from math import inf, log, floor

# --------------------------
# utilidades de grafo
# --------------------------
def add_edge(G, u, v, w):
    G[u].append((v, w))
    if v not in G:
        G[v] = G[v]  # garante v no dicionário

# --------------------------
# BaseCase (Algoritmo 2)
# S == {x} e x completo
# --------------------------
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
            if d_hat[u] + w <= d_hat[v] and d_hat[u] + w < B:  # ≤ conforme paper
                d_hat[v] = d_hat[u] + w
                if v not in inH:
                    heappush(H, (d_hat[v], v)); inH.add(v)
                else:
                    heappush(H, (d_hat[v], v))  # simplificação do DecreaseKey

    if len(U0) <= k:
        return B, U0
    else:
        Bp = max(d_hat[v] for v in U0)
        U = {v for v in U0 if d_hat[v] < Bp}
        return Bp, U

# --------------------------
# FindPivots (Algoritmo 1)
# --------------------------
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
            return set(S), W  # P = S

    # constrói floresta F de arestas em W que preservam d̂[v] = d̂[u]+w
    children = defaultdict(list)
    indeg = defaultdict(int)
    for u in W:
        for v, w in G[u]:
            if v in W and d_hat[u] + w == d_hat[v]:
                children[u].append(v)
                indeg[v] += 1

    # raízes em S com subárvore de tamanho >= k
    P = set()
    for u in S:
        if u not in W: 
            continue
        # BFS para contar (limitada por k)
        cnt = 0
        q = deque([u])
        seen = set([u])
        while q and cnt < k:
            x = q.popleft(); cnt += 1
            for y in children[x]:
                if y not in seen:
                    seen.add(y); q.append(y)
        if cnt >= k:
            P.add(u)
    return P, W

# --------------------------
# Estrutura D (Lema 3.3) - versão simples
# Mantém pares (key, value) e opera:
#   Insert, Pull (retira até M menores), BatchPrepend (insere valores menores que todos)
# --------------------------
class DStruct:
    def __init__(self, M, B):
        self.M = M
        self.B = B
        self.cur = {}           # key -> value
        self.small_buffer = []  # dados de BatchPrepend ainda não mesclados (lista simples)

    def insert(self, key, value):
        # mantém o menor valor por chave
        if key not in self.cur or value < self.cur[key]:
            self.cur[key] = value

    def batch_prepend(self, items):
        # itens são (key, value) todos < qualquer value atual (contrato do paper)
        # guardamos num buffer; na hora do pull, eles aparecem antes
        for k, v in items:
            if k not in self.cur or v < self.cur[k]:
                self.small_buffer.append((k, v))

    def pull(self):
        # materializa os menores até M
        # 1) junta buffer + cur em uma lista e pega M menores
        all_items = self.small_buffer + list(self.cur.items())
        if not all_items:
            return set(), self.B
        # ordena apenas o necessário (para didática)
        all_items.sort(key=lambda kv: kv[1])
        take = all_items[: self.M]
        S = {k for k, _ in take}
        # remove S de cur e do buffer
        cur_keys = set(self.cur.keys())
        for k in S:
            if k in self.cur:
                del self.cur[k]
        self.small_buffer = [(k, v) for (k, v) in self.small_buffer if k not in S]
        # calcula novo separador x (Bi do artigo)
        if self.small_buffer or self.cur:
            rest_min = min([v for (_, v) in self.small_buffer] + list(self.cur.values()))
            x = rest_min  # satisfaz max(S') < x ≤ min(D)
        else:
            x = self.B
        return S, x

# --------------------------
# BMSSP (Algoritmo 3) recursivo
# --------------------------
def bmssp(G, l, B, S, d_hat, t, k, K2_bound_factor=4):
    # l==0: base case
    if l == 0:
        return base_case(G, B, S, d_hat)

    P, W = find_pivots(G, B, S, d_hat)   # pivôs e conjunto W
    M = max(1, 2 ** ((l - 1) * t))       # M = 2^{(l-1)t}
    D = DStruct(M=M, B=B)
    for x in P:
        D.insert(x, d_hat[x])

    i = 0
    U = set()
    B0p = min([d_hat[x] for x in P], default=B)  # se P vazio
    last_Bp = B0p

    # limite de workload k*2^{l t} (usei 4k*2^{lt} para amortecer, como no Lema 3.9)
    U_cap = K2_bound_factor * (k * (2 ** (l * t)))

    while len(U) < U_cap:
        Si, Bi = D.pull()
        if not Si:
            # sucesso
            Bp = B
            Up = U | {x for x in W if d_hat[x] < Bp}
            return Bp, Up

        # chamada recursiva nível l-1
        Bp, Ui = bmssp(G, l - 1, Bi, Si, d_hat, t, k, K2_bound_factor)
        U |= Ui

        # relaxa arestas saindo de Ui e insere/antepõe em D conforme intervalos
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
        # devolve alguns de Si no intervalo [Bp, Bi)
        for x in Si:
            val = d_hat[x]
            if Bp <= val < Bi:
                Kbatch.append((x, val))
        if Kbatch:
            D.batch_prepend(Kbatch)

        last_Bp = Bp

    # encerramento parcial (workload grande)
    Up = U | {x for x in W if d_hat[x] < last_Bp}
    return last_Bp, Up

# --------------------------
# Função de alto nível: SSSP por BMSSP
# --------------------------
def sssp_break_sorting(G, source, B_limit=float("inf")):
    n = max(1, len(G))
    k = max(1, floor((log(n, 2)) ** (1/3)))               # ≈ log^{1/3} n
    t = max(1, floor((log(n, 2)) ** (2/3)))               # ≈ log^{2/3} n
    L = max(1, floor(log(n, 2) / max(1, t)))              # níveis

    # d̂ inicial
    d_hat = {u: inf for u in G}
    d_hat[source] = 0

    # S = {s}, l = L
    base_case.k = k
    find_pivots.k = k
    Bp, U = bmssp(G, L, B_limit, {source}, d_hat, t, k)

    # d̂ agora contém as distâncias para nós completados
    return d_hat, Bp, U, dict(meta=dict(n=n, k=k, t=t, L=L))


G = defaultdict(list)
add_edge(G,"A","B",2); add_edge(G,"A","C",5)
add_edge(G,"B","D",1); add_edge(G,"C","D",2)
add_edge(G,"D","E",3); add_edge(G,"C","F",10)
add_edge(G,"E","F",1)

dhat, Bp, U, meta = sssp_break_sorting(G, "A", B_limit=7)
print(meta, {v:dhat[v] for v in ["A","B","C","D","E","F"]})

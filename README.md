## BMSSP-Dijkstra

O BMSSP é uma alternativa eficiente ao Dijkstra para casos onde múltiplas fontes e um limite de distância são necessários. Ele mantém precisão local enquanto reduz custo computacional, sendo ideal para aplicações em redes grandes como mapas, redes sociais e infraestruturas de transporte.

## Areas de atuação
BMSSP é eficiente quando múltiplas fontes e um limite de distância são necessários.
Mantém precisão local, reduzindo custo computacional.
Aplicações: mapas, redes sociais, transporte.

Superando o algoritmo de Dijkstra em diversos casos

Fonte do artigo : https://arxiv.org/pdf/2504.17033 
 
BMSSP×Dijkstra - Comparativo entre 2 algoritmos de rotas

## GRAPOS.PY
python grafos.py

# Retorno
{'meta': {'n': 6, 'k': 1, 't': 1, 'L': 2}} {'A': 0, 'B': 2, 'C': 5, 'D': 3, 'E': 6, 'F': 7}

# GRAPOS2.PY
python grafos2.py

# Retorno
python3 grafos2.py
== BMSSP vs Dijkstra ==
Tempo BMSSP:   0.000601 s
Tempo Dijkstra:1.2e-05 s
Vértices verificados (dist <= 7): 6
Inconsistências: 0
Metadados BMSSP: {'n': 6, 'k': 1, 't': 1, 'L': 2, 'Bp': 7, 'U_size': 5, 'time_s': 0.0006013329839333892}
A: BMSSP=0 | Dijkstra=0
B: BMSSP=2 | Dijkstra=2
C: BMSSP=5 | Dijkstra=5
D: BMSSP=3 | Dijkstra=3
E: BMSSP=6 | Dijkstra=6
F: BMSSP=7 | Dijkstra=7

# Descrição do Dataset SNAP roadNet-CA
Esse dataset faz parte da coleção SNAP (Stanford Network Analysis Platform). Veja os detalhes:

É um grafo não-direcionado que modela a rede viária da Califórnia.
Os nós representam interseções e pontos finais de estradas, enquanto as arestas conectam esses pontos
Estatísticas principais:
Nº de nós: 1 965 206
Nº de arestas: 5 533 214 (ou 2 766 607 no SNAP, considerado como sem duplicação)

A maior componente conectada (WCC) inclui ~99,6% dos nós (~1 957 027) e ~99,8% das arestas (~5 520 776)
cise.ufl.edu
snap.stanford.edu

Coeficiente médio de clustering: ~0,0464 (é baixo)
cise.ufl.edu
snap.stanford.edu
Número de triângulos fechados: ~120 676
Fração de triângulos fechados: ~0,02–0,06 (varia por fonte)

sparse.tamu.edu
cise.ufl.edu
snap.stanford.edu
studentwork.prattsi.org
deepai.org
researchgate.net

Diâmetro (caminho mais longo entre dois nós): aproximadamente 849–850

Diâmetro efetivo (percentil 90%): ~500
snap.stanford.edu
cise.ufl.edu

Contexto:
Publicado em 2008 por Jure Leskovec, K. Lang, A. Dasgupta e M. Mahoney, parte do estudo sobre estruturas comunitárias em grandes redes disponível em Internet Mathematics.
networks.skewed.de
snap.stanford.edu
sparse.tamu.edu

Disponível em diversos formatos (texto .txt.gz, Matrix Market, SuiteSparse, etc.) e comumente usado para benchmarking e experimentos com algoritmos de grafos e roteamento.

# Versão que utliza dataset com rotas da Florida bmssp_snap_dimacs

pip install requests tqdm

python bmssp_snap_dimacs.py

Por padrão roda a demo SNAP (roadNet-CA), baixa o .txt.gz, constrói o grafo (peso=1) e compara BMSSP×Dijkstra com validação por amostra (configurável).

Para testar DIMACS Florida (ponderado), abra o arquivo e descomente a linha demo_dimacs(...).

Parâmetros principais dentro do script
B_limit: raio/limite de custo para o BMSSP (ex.: 800 em SNAP; 200000 em DIMACS).

sample_check: tamanho da amostra para validação contra Dijkstra em grafos grandes (reduz custo de checagem).

seed: controla a escolha da fonte (source) aleatória.

# Retorno 
== Relatório SNAP ==
{'bmssp_time_s': 1.883903, 
'dijkstra_time_s': 6.485129, 
'checked_vertices': 19930, 
'mismatches': 19906, 
'meta': {'n': 1965206, 'k': 2, 't': 7, 'L': 2, 'Bp': 800, 'U_size': 45517, 'time_s': 1.8839032619725913}, 'examples': [(132327, inf, 399), (552897, inf, 160), (176340, inf, 466), (1611459, inf, 403), (867363, inf, 312), (550408, inf, 178), (216808, inf, 354), (1758830, inf, 382), (1893278, inf, 345), (1834731, inf, 549)]}


# BMSSP-Dijkstra
BMSSP×Dijkstra - Comparativo entre 2 algoritmos de rotas


pip install requests tqdm


python bmssp_snap_dimacs.py

Por padrão roda a demo SNAP (roadNet-CA), baixa o .txt.gz, constrói o grafo (peso=1) e compara BMSSP×Dijkstra com validação por amostra (configurável).

Para testar DIMACS Florida (ponderado), abra o arquivo e descomente a linha demo_dimacs(...).

Parâmetros principais dentro do script
B_limit: raio/limite de custo para o BMSSP (ex.: 800 em SNAP; 200000 em DIMACS).

sample_check: tamanho da amostra para validação contra Dijkstra em grafos grandes (reduz custo de checagem).

seed: controla a escolha da fonte (source) aleatória.

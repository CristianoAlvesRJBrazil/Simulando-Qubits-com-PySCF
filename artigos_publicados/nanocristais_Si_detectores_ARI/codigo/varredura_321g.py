#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""varredura_321g.py -- teste de robustez da sensibilidade ao criterio de corte:
repete, na base de valencia dividida 3-21G, as comparacoes entre criterios que a
Tabela 2 traz em STO-3G. Responde a pergunta "o efeito do criterio sobrevive a
uma base adequada?". Executa na GPU (SCF direto, gpu4pyscf); grava em ../dados/gaps.json.
Ja disponiveis em 3-21G (nao recalcula): Si10H16 (0,8 nm, mv2) e Si22H28 (1,0 nm, mv2).
Uso: python3 varredura_321g.py"""
import time
import varredura_gap as vg
from varredura_extra import construir_preencher
from varredura_gpu import calcular_gpu

# do mais barato ao mais caro (N_b = 13 n_Si + 2 n_H)
SISTEMAS = [
    (0.8, 'mv1',  'min_viz<=1',          lambda: vg.construir(0.8, 8, 1)),   # Si11H18, 179
    (1.0, 'mv1',  'min_viz<=1',          lambda: vg.construir(1.0, 8, 1)),   # Si26H36, 410
    (1.2, 'mv2',  'min_viz=2',           lambda: vg.construir(1.2, 8, 2)),   # Si36H40, 548
    (1.2, 'mv2p', 'min_viz=2+preencher', lambda: construir_preencher(1.2)),  # Si39H40, 587
    (1.2, 'mv1',  'min_viz<=1',          lambda: vg.construir(1.2, 8, 1)),   # Si45H58, 701
]

if __name__ == '__main__':
    db = vg.carregar(); t0 = time.time()
    for D, tag, crit, gerar in SISTEMAS:
        si, h = gerar()
        print("== %s  D=%.1f nm  %s  (N_b=%d) ==" % (vg.formula(si, h), D, crit, 13*len(si) + 2*len(h)), flush=True)
        for met in ('HF', 'PBE'):
            calcular_gpu(db, 'b321', D, tag, crit, si, h, '3-21g', met, direto=True)
        vg.exportar_csv(db)
    print("FIM. Total %.0f s" % (time.time() - t0), flush=True)

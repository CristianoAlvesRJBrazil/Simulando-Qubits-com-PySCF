#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relaxacao de CONTROLE em base 3-21G (PBE), reutilizando relaxacao.py.
Objetivo: separar o efeito fisico da relaxacao do artefato de contracao de ligacoes da base minima.
Grava em ../dados/relaxacao_321g.json e XYZ com prefixo '321g_'. Uso: python3 relaxacao_321g.py Si10H16 [Si19H28]"""
import os, sys, time
import relaxacao as rx
rx.BASE = "3-21g"
rx.JSON_PATH = os.path.join(rx.DADOS, "relaxacao_321g.json")
_gravar = rx.gravar_xyz
def gravar_321g(nome, si, h, comentario=""):
    _gravar("321g_" + nome.replace("sto3g", "321g"), si, h, comentario.replace("STO-3G", "3-21G"))
rx.gravar_xyz = gravar_321g
if __name__ == "__main__":
    rotulos = sys.argv[1:] or ["Si10H16"]
    for r in rotulos:
        t0 = time.time()
        reg = rx.rodar(r)
        reg["tempo_total_s"] = time.time() - t0
        dados = rx.carregar_json(); dados["_meta"]["base"] = "3-21g"
        dados["sistemas"][r] = reg
        rx.salvar_json(dados)
        print("  [%s gravado em %s; %.0f s]" % (r, rx.JSON_PATH, reg["tempo_total_s"]))

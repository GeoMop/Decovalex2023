# Decovalex2023

Tool for upsaling of a Discrete Fracture Network and repository sources into eqivalent 
heterogeneous continuum on a structed PFlotran grid.

Based on original SANDIA scripts (see below).


## Instalation
Installation to user space (no root / administrator rights required). 

From local sources in a `decovalex` directory:

    ```
    python3 -m pip install ./decovalex
    ```
You can clone the git repository (get sources) by:
    
    ```
    git clone https://github.com/GeoMop/Decovalex2023.git
    ```

Directly from the github repository:
    ```
    python3 -m pip install "decovalex @ git+https://github.com/GeoMop/Decovalex2023
    ``` 
    


## Repository structure
- `decodfn` 
    - `main.py` : the main script, the `main()` called by the command `decodfn`
- `tests` : test calculations 
    - `dfn_253` : test sample of 253 fractures
    - `repo_multiscale` : basic multiscale repository model
    
   
1) skript na promítnutí DFN do heterogenní propustnosti v pravidelné mřížce (porous medium - CPM)
mapDFN-krychlicka.zip
Zde je to v přehlednější verzi se závislostmi a s celými vstupními a výstupními daty, je to na jednodušší geometrii krychličky s puklinami, hrubá diskretizace s malou velikostí souborů (není to ale přímo ke spuštění protože to není ve správných cestách). Hlavní skript je mapdfn2pflotran.py, který obsahuje přímo v kódu veškeré vstupy na prvních 10-20 řádcích.
Vstupem jsou .txt soubory s daty puklin generovaných DFNworks (to už tam komplet nedávám, teď nebylo relevantní).
Výstupem jsou .h5 soubory.

uloziste-mapDFN+pflotran_krok25.zip
Zde je to pro aktuální úlohu celého úložiště. Nemám kompletní vstupní data. Skript je trochu upravený, proto jsem tím byl zmatený když jsem to ukazoval - teď jsem při podrobném čtení zjistil, že část dat se dává jako argumenty příkazové řádky. Zároveň je tam ale x různých specifických nastavení pro tu konkrétní úlohu, např. varianta “repo” se zjemněním. Podle komentářů to vypadá, že to aktuálně dělal někdo z DoE pro potřeby této úlohy D2023.

2) úloha Pflotran
Vstup pro Pflotran je soubor .in a všechny generované .h5 z mapdfn (kromě řídícího mapELLIPSES viz níže). Pro aktuální úlohu je součástí zipu výše. Pro úlohu krychlička je tento samostatný Pflotran-krychlicka.zip
V tom .in pro celé úložiště je definovaný region “dgr” a v něm je zadaná okrajová podmínka “waste form” což je speciální modul na zdrojový člen (samotný kontejner). Teď je to ted udělané tak, že jakoby byl obří kontejner přes cca 100m blok buněk rovnou v žule (to je stav do kterého se dostali než se obrátili na mě)..

3) zobrazení v ParaView
Soubory .H5 mají mezi sebou závislosti (na to jsem zapomněl). Je hlavní soubor mapELLIPSES.h5 a ten se odkazuje na zbylé:
anisotropic_k.h5
isotropic_k.h5
materials.h5
porosity.h5
tortuosity.h5
V paraview je třeba otevřít ten mapellipses. Načte si nějak ty zbylé a v menu načtených dat to rovnou nabídne k zaškrtnutí pole veličin, které se mají načíst (K, poro, index pukliny,...). Pak už se běžně ještě přepne zobrazení že má být surface a že colouring pole toho datového pole. 

4) Adresář se skriptem pro generování .h5 pro okrajovou podmínku. Jak je to přesně dělané jsem nestudoval. To si vytvářel sám Ondrej Mikláš.

## Original scripts authors

Authors: Emily Stein (ergiamb@sandia.gov)
        Kris Kuhlman (klkuhl@sandia.gov)
        Applied Systems Analysis and Research, 8844
        Sandia National Laboratories

Date: 07/13/18
SAND Number: SAND2018-7604 O

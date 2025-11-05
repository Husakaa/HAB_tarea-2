#!/usr/bin/env python3
"""
Network Propagation - GUILD & DIAMOnD
======================================

Propagación en redes de interacción proteína-proteína para identificar
genes relacionados con genes semilla.

ENTRADAS:
    -n: Archivo de red
        • data/network_guild.txt (formato GUILD)
        • data/network_diamond.txt (formato DIAMOnD) 
        • data/string_network_filtered_hugo-400.tsv (formato STRING)
    
    -g: Genes semilla (manual): ENO1 PGK1 HK2
    -s: Archivo de genes: data/genes_seed.txt

SALIDAS:
    results/guild_results.tsv
    results/diamond_results.tsv

EJEMPLOS:
    python network_propagation.py -n data/network_guild.txt -g 2023 5230 3099
    python network_propagation.py -n data/string_network_filtered_hugo-400.tsv -s data/genes_seed.txt
"""

import argparse
import os
import sys
import networkx as nx
import numpy as np
import pandas as pd
import mygene
from scipy.stats import hypergeom


def cargar_red(archivo):
    """Carga red detectando formato automáticamente"""
    with open(archivo, 'r') as f:
        primera = f.readline().strip()
        
    if '\t' in primera or 'protein' in primera.lower():
        # STRING: TSV
        print("   Formato: STRING (TSV)")
        df = pd.read_csv(archivo, sep='\t')
        return nx.from_pandas_edgelist(df, df.columns[0], df.columns[1])
    
    elif ',' in primera:
        # DIAMOnD: CSV (puede estar en .txt o .csv)
        print("   Formato: DIAMOnD (CSV)")
        df = pd.read_csv(archivo, header=None, names=['source', 'target'])
        return nx.from_pandas_edgelist(df, 'source', 'target')
    
    else:
        # GUILD: nodo1 peso nodo2 (espacios)
        partes = primera.split()
        
        if len(partes) == 3:
            print("   Formato: GUILD")
            G = nx.Graph()
            with open(archivo, 'r') as f:
                for linea in f:
                    partes = linea.strip().split()
                    if len(partes) == 3:
                        G.add_edge(partes[0], partes[2], weight=float(partes[1]))
            return G
        else:
            raise ValueError(f"Formato no reconocido. Primera línea: {primera}")


def leer_genes(archivo):
    """Lee genes desde archivo"""
    with open(archivo, 'r') as f:
        contenido = f.read().replace(' y ', ' ').replace(',', ' ')
    return [g.strip() for g in contenido.split() if g.strip()]


def agregar_genes_a_red(archivo_red, genes_faltantes, archivo_string=None):
    """
    Agrega genes faltantes a la red buscando sus interacciones en STRING.
    
    Args:
        archivo_red: Archivo de red GUILD/DIAMOnD a modificar
        genes_faltantes: Lista de Entrez IDs que faltan
        archivo_string: Archivo STRING de donde obtener interacciones (opcional)
    
    Returns:
        bool: True si se agregaron genes, False si no
    """
    # Buscar archivo STRING en múltiples ubicaciones
    posibles_rutas = [
        archivo_string,
        'data/string_network_filtered_hugo-400.tsv',
    ]
    
    archivo_string_encontrado = None
    for ruta in posibles_rutas:
        if ruta and os.path.exists(ruta):
            archivo_string_encontrado = ruta
            break
    
    if not archivo_string_encontrado:
        print("   Advertencia: No se encontró archivo STRING, no se pueden agregar genes")
        return False
    
    # Convertir Entrez IDs a símbolos usando mygene
    try:
        mg = mygene.MyGeneInfo()
        
        print("   Convirtiendo Entrez IDs a símbolos...")
        results = mg.querymany(genes_faltantes, scopes='entrezgene', fields='symbol', species='human')
        
        # Crear mapeo de Entrez -> símbolo
        entrez_a_simbolo = {}
        simbolo_a_entrez = {}
        
        for r in results:
            entrez = str(r.get('query'))
            simbolo = r.get('symbol')
            if simbolo:
                entrez_a_simbolo[entrez] = simbolo
                simbolo_a_entrez[simbolo] = entrez
        
        if not entrez_a_simbolo:
            print("   No se pudieron convertir los Entrez IDs")
            return False
        
        simbolos_buscar = list(entrez_a_simbolo.values())
        
    except ImportError:
        print("   Error: mygene no instalado, no se pueden convertir IDs")
        return False
    except Exception as e:
        print(f"   Error al convertir IDs: {e}")
        return False
    
    print(f"   Buscando interacciones en STRING para: {', '.join(simbolos_buscar)}...")
    
    # Leer red STRING
    try:
        df_string = pd.read_csv(archivo_string_encontrado, sep='\t')
    except Exception as e:
        print(f"   Error al leer STRING: {e}")
        return False
    
    # Buscar interacciones que involucren estos genes
    interacciones = df_string[
        (df_string['protein1_hugo'].isin(simbolos_buscar)) | 
        (df_string['protein2_hugo'].isin(simbolos_buscar))
    ]
    
    # Filtrar solo interacciones entre los genes de interés
    interacciones_internas = interacciones[
        (interacciones['protein1_hugo'].isin(simbolos_buscar)) & 
        (interacciones['protein2_hugo'].isin(simbolos_buscar))
    ]
    
    if len(interacciones_internas) == 0:
        print("   No se encontraron interacciones entre estos genes en STRING")
        return False
    
    print(f"   Encontradas {len(interacciones_internas)} interacciones")
    
    # Preparar líneas a añadir
    nuevas_lineas = []
    interacciones_unicas = set()
    
    for _, row in interacciones_internas.iterrows():
        gene1 = row['protein1_hugo']
        gene2 = row['protein2_hugo']
        
        if gene1 in simbolo_a_entrez and gene2 in simbolo_a_entrez:
            entrez1 = simbolo_a_entrez[gene1]
            entrez2 = simbolo_a_entrez[gene2]
            
            # Evitar duplicados (orden no importa)
            par = tuple(sorted([entrez1, entrez2]))
            if par not in interacciones_unicas:
                interacciones_unicas.add(par)
                
                # Formato según tipo de red
                if 'guild' in archivo_red.lower():
                    nuevas_lineas.append(f"{entrez1} 1 {entrez2}")
                else:  # diamond
                    nuevas_lineas.append(f"{entrez1},{entrez2}")
    
    if not nuevas_lineas:
        return False
    
    # Añadir líneas al archivo
    try:
        with open(archivo_red, 'a') as f:
            for linea in nuevas_lineas:
                f.write('\n' + linea)
        
        print(f"   Agregadas {len(nuevas_lineas)} interacciones a {archivo_red}")
        return True
    except Exception as e:
        print(f"   Error al escribir en {archivo_red}: {e}")
        return False


def convertir_a_entrez(genes):
    """
    Convierte símbolos de genes a Entrez IDs usando MyGene.info
    
    Args:
        genes: Lista de símbolos de genes (ej: ['ENO1', 'PGK1', 'HK2'])
    
    Returns:
        dict: Mapeo {símbolo: entrez_id}
    """
    print("\nConvirtiendo símbolos a Entrez IDs...")
    
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        
        # Consultar MyGene
        results = mg.querymany(genes, scopes='symbol', fields='entrezgene', species='human')
        
        # Crear mapeo
        mapeo = {}
        no_encontrados = []
        
        for r in results:
            simbolo = r.get('query')
            entrez = r.get('entrezgene')
            
            if entrez:
                mapeo[simbolo] = str(int(entrez))
                print(f"   {simbolo} -> {int(entrez)}")
            else:
                no_encontrados.append(simbolo)
        
        if no_encontrados:
            print(f"\nAdvertencia: No se encontró Entrez ID para: {', '.join(no_encontrados)}")
        
        return mapeo
        
    except ImportError:
        print("\nError: mygene no instalado.")
        print("Instala con: pip install mygene")
        sys.exit(1)
    except Exception as e:
        print(f"\nError al consultar MyGene.info: {e}")
        print("Verifica tu conexión a internet")
        sys.exit(1)


def guild(G, semillas):
    """
    GUILD usando PageRank personalizado de NetworkX
    (equivalente a Random Walk with Restart)
    """
    print("\nEjecutando GUILD (PageRank)...")
    
    # Personalización: probabilidad 1 en semillas, 0 en el resto
    personalizacion = {nodo: 1.0 if nodo in semillas else 0.0 for nodo in G.nodes()}
    
    # PageRank personalizado
    scores = nx.pagerank(G, personalization=personalizacion, alpha=0.85)
    
    # Crear DataFrame
    df = pd.DataFrame([
        {'gene': nodo, 'guild_score': score, 'is_seed': nodo in semillas}
        for nodo, score in scores.items()
    ]).sort_values('guild_score', ascending=False)
    
    return df


def diamond(G, semillas, max_genes=200):
    """
    DIAMOnD: Disease Module Detection
    Algoritmo iterativo con test hipergeométrico
    """
    print("\nEjecutando DIAMOnD...")
    
    modulo = set(semillas)
    candidatos = set(G.nodes()) - modulo
    n_total = G.number_of_nodes()
    resultados = []
    
    for i in range(min(max_genes, len(candidatos))):
        mejor = None
        mejor_pval = 1.0
        
        # Evaluar todos los candidatos
        for candidato in candidatos:
            conexiones = sum(1 for vecino in G.neighbors(candidato) if vecino in modulo)
            grado = G.degree(candidato)
            
            if grado > 0:
                # Test hipergeométrico
                pval = hypergeom.sf(conexiones - 1, n_total, len(modulo), grado)
                
                if pval < mejor_pval:
                    mejor_pval = pval
                    mejor = (candidato, conexiones, grado)
        
        if not mejor:
            break
        
        # Añadir mejor candidato al módulo
        gen, conex, grado = mejor
        modulo.add(gen)
        candidatos.remove(gen)
        
        resultados.append({
            'rank': i + 1,
            'gene': gen,
            'p_value': mejor_pval,
            'diamond_score': -np.log10(mejor_pval + 1e-300),
            'connections_to_module': conex,
            'total_degree': grado
        })
        
        if (i + 1) % 50 == 0:
            print(f"   Procesados {i + 1}/{max_genes} genes...")
    
    df = pd.DataFrame(resultados)
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Propagación en redes - GUILD & DIAMOnD',
        epilog="""
Ejemplos:
  # Uso básico (usa redes por defecto)
  python %(prog)s -g 2023 5230 3099
  
  # Con red personalizada
  python %(prog)s -n data/string_network_filtered_hugo-400.tsv -s data/genes_seed.txt
  
  # Solo un algoritmo
  python %(prog)s -g 2023 5230 3099 -a guild
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-n', '--network',
                       help='Archivo de red (opcional, usa data/network_guild.txt y data/network_diamond.txt por defecto)')
    parser.add_argument('-g', '--genes', nargs='+', 
                       help='Genes semilla (ej: 2023 5230 3099 para Entrez IDs)')
    parser.add_argument('-s', '--seed-file', 
                       help='Archivo con genes semilla')
    parser.add_argument('-a', '--algorithm', 
                       choices=['guild', 'diamond', 'both'], 
                       default='both',
                       help='Algoritmo a ejecutar (default: both)')
    parser.add_argument('-o', '--output', 
                       default='results',
                       help='Directorio de salida (default: results)')
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("  PROPAGACIÓN EN REDES - GUILD & DIAMOnD")
    print("="*70)
    
    # Determinar genes semilla
    if args.seed_file:
        if not os.path.exists(args.seed_file):
            print(f"\nError: {args.seed_file} no encontrado")
            sys.exit(1)
        genes = leer_genes(args.seed_file)
        print(f"\nGenes desde {args.seed_file}: {', '.join(genes)}")
    
    elif args.genes:
        genes = args.genes
        print(f"\nGenes semilla: {', '.join(genes)}")
    
    else:
        # Buscar automáticamente
        for path in ['genes_seed.txt', 'data/genes_seed.txt']:
            if os.path.exists(path):
                genes = leer_genes(path)
                print(f"\nGenes desde {path}: {', '.join(genes)}")
                break
        else:
            genes = ['2023', '5230', '3099']  # Entrez IDs por defecto
            print(f"\nGenes por defecto: {', '.join(genes)} (Entrez IDs)")
    
    # Crear directorio de salida
    os.makedirs(args.output, exist_ok=True)
    
    # Si no se especifica red, usar redes por defecto según algoritmo
    if not args.network:
        print("\nUsando redes por defecto (requieren Entrez IDs):")
        
        # Detectar si los genes son símbolos y convertir
        genes_convertidos = genes
        es_simbolo = any(not g.isdigit() for g in genes)
        
        if es_simbolo:
            print("Genes proporcionados parecen símbolos, convirtiendo a Entrez IDs...")
            mapeo = convertir_a_entrez(genes)
            genes_convertidos = [mapeo[g] for g in genes if g in mapeo]
            
            if not genes_convertidos:
                print("Error: No se pudieron convertir los genes")
                sys.exit(1)
        
        # Aplicar GUILD
        if args.algorithm in ['guild', 'both']:
            network_guild = 'data/network_guild.txt'
            if not os.path.exists(network_guild):
                print(f"Error: {network_guild} no encontrado")
                sys.exit(1)
            
            print(f"\nGUILD - Cargando red: {network_guild}")
            G_guild = cargar_red(network_guild)
            print(f"   {G_guild.number_of_nodes()} nodos, {G_guild.number_of_edges()} aristas")
            
            validos_guild = [g for g in genes_convertidos if g in G_guild.nodes()]
            if not validos_guild:
                print("Error: Ningún gen en la red GUILD")
                sys.exit(1)
            print(f"Genes validos: {', '.join(validos_guild)}")
            
            resultados_guild = guild(G_guild, validos_guild)
            archivo_guild = f"{args.output}/guild_results.tsv"
            resultados_guild.to_csv(archivo_guild, sep='\t', index=False)
            print(f"Guardado: {archivo_guild}")
        
        # Aplicar DIAMOND
        if args.algorithm in ['diamond', 'both']:
            network_diamond = 'data/network_diamond.txt'
            if not os.path.exists(network_diamond):
                print(f"Error: {network_diamond} no encontrado")
                sys.exit(1)
            
            print(f"\nDIAMOnD - Cargando red: {network_diamond}")
            G_diamond = cargar_red(network_diamond)
            print(f"   {G_diamond.number_of_nodes()} nodos, {G_diamond.number_of_edges()} aristas")
            
            validos_diamond = [g for g in genes_convertidos if g in G_diamond.nodes()]
            if not validos_diamond:
                print("Error: Ningún gen en la red DIAMOnD")
                sys.exit(1)
            print(f"Genes validos: {', '.join(validos_diamond)}")
            
            resultados_diamond = diamond(G_diamond, validos_diamond)
            archivo_diamond = f"{args.output}/diamond_results.tsv"
            resultados_diamond.to_csv(archivo_diamond, sep='\t', index=False)
            print(f"Guardado: {archivo_diamond}")
    
    else:
        # Red personalizada especificada con -n
        print(f"\nCargando red: {args.network}")
        if not os.path.exists(args.network):
            print("Error: Archivo no encontrado")
            sys.exit(1)
        
        G = cargar_red(args.network)
        print(f"   {G.number_of_nodes()} nodos, {G.number_of_edges()} aristas")
        
        # Verificar genes en la red
        validos = [g for g in genes if g in G.nodes()]
        no_encontrados = [g for g in genes if g not in G.nodes()]
        
        # Si hay genes no encontrados y parecen símbolos, intentar conversión
        if no_encontrados and any(not g.isdigit() for g in no_encontrados):
            # Detectar si es red GUILD/DIAMOND (requiere Entrez IDs)
            if 'guild' in args.network.lower() or 'diamond' in args.network.lower():
                print(f"\nGenes no encontrados: {', '.join(no_encontrados)}")
                print("Intentando conversión automática a Entrez IDs...")
                
                mapeo = convertir_a_entrez(no_encontrados)
                
                # Verificar si los IDs convertidos están en la red
                ids_convertidos_faltantes = []
                genes_finales = []
                
                for gen in genes:
                    if gen in mapeo:
                        entrez_id = mapeo[gen]
                        if entrez_id in G.nodes():
                            genes_finales.append(entrez_id)
                            print(f"   {gen} ({entrez_id}) encontrado en red")
                        else:
                            ids_convertidos_faltantes.append(entrez_id)
                            print(f"   {gen} ({entrez_id}) NO encontrado en red")
                    elif gen in G.nodes():
                        genes_finales.append(gen)
                
                # Si hay IDs convertidos que no están en la red, intentar agregarlos
                if ids_convertidos_faltantes:
                    print("\nIntentando agregar genes faltantes a la red...")
                    if agregar_genes_a_red(args.network, ids_convertidos_faltantes):
                        # Recargar red
                        print("\nRecargando red actualizada...")
                        G = cargar_red(args.network)
                        print(f"   {G.number_of_nodes()} nodos, {G.number_of_edges()} aristas")
                        
                        # Volver a verificar
                        genes_finales = []
                        for gen in genes:
                            if gen in mapeo and mapeo[gen] in G.nodes():
                                genes_finales.append(mapeo[gen])
                            elif gen in G.nodes():
                                genes_finales.append(gen)
                
                validos = genes_finales
                no_encontrados = [g for g in genes if g not in mapeo or mapeo.get(g) not in G.nodes()]
        
        if no_encontrados and len(validos) == 0:
            print(f"\nGenes no encontrados en red (después de conversión): {', '.join(no_encontrados)}")
        
        if not validos:
            print("\nError: Ningún gen semilla en la red")
            sys.exit(1)
        
        print(f"\nGenes validos: {', '.join(validos)}")
        
        # Ejecutar algoritmos
        if args.algorithm in ['guild', 'both']:
            resultados_guild = guild(G, validos)
            archivo_guild = f"{args.output}/guild_results.tsv"
            resultados_guild.to_csv(archivo_guild, sep='\t', index=False)
            print(f"\nGuardado: {archivo_guild}")
        
        if args.algorithm in ['diamond', 'both']:
            resultados_diamond = diamond(G, validos)
            archivo_diamond = f"{args.output}/diamond_results.tsv"
            resultados_diamond.to_csv(archivo_diamond, sep='\t', index=False)
            print(f"\nGuardado: {archivo_diamond}")
    
    print("\n" + "="*70)
    print("ANALISIS COMPLETADO")
    print("="*70)
    print(f"\nResultados en: {args.output}/\n")


if __name__ == '__main__':
    main()
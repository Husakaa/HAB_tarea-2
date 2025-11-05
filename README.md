# Propagación en Redes de Interacción Proteína-Proteína

Pipeline completo en Python para identificar genes candidatos relacionados con genes semilla utilizando algoritmos de propagación en redes: **GUILD** y **DIAMOnD**.

## Descripción

Este proyecto implementa dos algoritmos complementarios para el análisis de redes de interacción proteína-proteína (PPI):

- **GUILD** (Genes Underlying Inheritance Linked Disorders) - Propagación mediante Random Walk with Restart
- **DIAMOnD** (DIseAse Module Detection) - Detección de módulos de enfermedad mediante significancia estadística

El análisis permite identificar genes candidatos que están funcionalmente relacionados con un conjunto de genes semilla conocidos.

## Instalación (recomendado con entorno virtual)

1. `python -m venv .venv`
2. `source .venv/bin/activate`  (Linux/Mac) ó 
   `.venv\Scripts\activate`     (Windows)
3. `pip install -r requirements.txt`

### Requisitos

- Python 3.7+
- R 4.0+
- RStudio (recomendado)
 
### Dependencias Python
```bash
pip install -r requirements.txt
```

## Estructura del Proyecto
```
proyecto/
├── data/
│   ├── genes_seed.txt                          # Genes semilla (entrada)
│   ├── network_guild.txt                       # Red PPI formato GUILD (Entrez IDs)
│   ├── network_diamond.txt                     # Red PPI formato DIAMOnD (Entrez IDs)
│   └── string_network_filtered_hugo-400.tsv    # Red STRING (símbolos HUGO)
├── scripts/
│   └── network_propagation.py                  # Script principal
├── results/
│   ├── guild_results.tsv                       # Resultados GUILD
│   └── diamond_results.tsv                     # Resultados DIAMOnD
├── requirements.txt
└── README.md
```

## Formatos de Red Soportados

El script detecta automáticamente el formato de red:

### 1. **GUILD** (formato: `nodo1 peso nodo2`)
```
2023 1 5230
5230 1 3099
2023 0.8 3099
```
- Usa **Entrez IDs** como identificadores
- Pesos opcionales (default: 1)
- Separado por espacios

### 2. **DIAMOnD** (formato: `nodo1,nodo2`)
```
2023,5230
5230,3099
2023,3099
```
- Usa **Entrez IDs** como identificadores
- Separado por comas
- Puede estar en archivo `.txt` o `.csv`

### 3. **STRING** (formato: TSV con encabezados)
```
protein1_hugo    protein2_hugo    combined_score
ENO1            PGK1             950
PGK1            HK2              850
```
- Usa **símbolos HUGO** como identificadores
- Formato tabular (TSV)
- Incluye scores de confianza

## Uso

### Ejecución Básica

```bash
# Usar redes por defecto (GUILD y DIAMOnD) con genes desde archivo
python scripts/network_propagation.py -s data/genes_seed.txt

# Genes directos (símbolos HUGO)
python scripts/network_propagation.py -g ENO1 PGK1 HK2

# Genes directos (Entrez IDs)
python scripts/network_propagation.py -g 2023 5230 3099
```

### Opciones Avanzadas

```bash
# Red personalizada (STRING con símbolos)
python scripts/network_propagation.py -n data/string_network_filtered_hugo-400.tsv -s data/genes_seed.txt

# Red personalizada (GUILD/DIAMOnD con Entrez IDs)
python scripts/network_propagation.py -n data/network_guild.txt -g 2023 5230 3099

# Ejecutar solo GUILD
python scripts/network_propagation.py -g ENO1 PGK1 HK2 -a guild

# Ejecutar solo DIAMOnD
python scripts/network_propagation.py -g ENO1 PGK1 HK2 -a diamond

# Especificar directorio de salida
python scripts/network_propagation.py -s data/genes_seed.txt -o mis_resultados
```

### Argumentos Disponibles

| Argumento | Descripción | Ejemplo |
|-----------|-------------|---------|
| `-g, --genes` | Genes semilla (separados por espacio) | `-g ENO1 PGK1 HK2` |
| `-s, --seed-file` | Archivo con genes semilla | `-s data/genes_seed.txt` |
| `-n, --network` | Archivo de red personalizado | `-n data/network_guild.txt` |
| `-a, --algorithm` | Algoritmo a ejecutar (`guild`/`diamond`/`both`) | `-a guild` |
| `-o, --output` | Directorio de salida | `-o results` |

## Formato de Genes de Entrada

### Opción 1: Archivo de texto (`genes_seed.txt`)

**Formato flexible:**
```
ENO1, PGK1 y HK2
```

O un gen por línea:
```
ENO1
PGK1
HK2
```

### Opción 2: Línea de comandos
```bash
python scripts/network_propagation.py -g ENO1 PGK1 HK2
```

## Conversión Automática de IDs

El script convierte automáticamente entre **símbolos HUGO** y **Entrez IDs** usando MyGene:

```
ENO1 → 2023
PGK1 → 5230
HK2  → 3099
```

**Reglas:**
- Redes GUILD/DIAMOnD → requieren **Entrez IDs**
- Red STRING → usa **símbolos HUGO**
- Conversión automática según el tipo de red

## Resultados

### Archivos Generados

#### `guild_results.tsv`
```
rank    gene    score
1       ENO1    0.168017
2       PGK1    0.167927
3       HK2     0.167544
4       GAPDH   0.001569
5       TPI1    0.001419
```

Columnas:
- **rank**: Posición del gen (1 = más relevante)
- **gene**: Identificador del gen
- **score**: Score de propagación (mayor = más cercano a genes semilla)

#### `diamond_results.tsv`
```
rank    gene    p_value
1       GAPDH   0.000123
2       TPI1    0.000456
3       ALDOA   0.001234
```

Columnas:
- **rank**: Posición del gen
- **gene**: Identificador del gen
- **p_value**: Significancia estadística (menor = más significativo)

### Interpretación de Resultados

#### GUILD
- **Score alto (>0.1)**: Genes semilla originales
- **Score medio (0.001-0.1)**: Genes muy conectados con semillas
- **Score bajo (<0.001)**: Genes distantes en la red

**Recomendación**: Usar genes con score > 0.001

#### DIAMOnD
- **p < 0.001**: Alta significancia (muy recomendado)
- **p < 0.05**: Significativo (recomendado)
- **p > 0.05**: No significativo (descartar)

**Recomendación**: Filtrar resultados con p-value < 0.05

## Ejemplos Completos

### Ejemplo 1: Análisis de genes glucolíticos con STRING
```bash
# 1. Crear archivo de genes
echo "ENO1, PGK1, HK2" > data/genes_seed.txt

# 2. Ejecutar análisis con red STRING
python scripts/network_propagation.py \
  -n data/string_network_filtered_hugo-400.tsv \
  -s data/genes_seed.txt

# 3. Ver resultados
cat results/guild_results.tsv | head -20
cat results/diamond_results.tsv | head -20
```

### Ejemplo 2: Análisis con redes por defecto
```bash
# Usar Entrez IDs directamente
python scripts/network_propagation.py -g 2023 5230 3099

# Resultados en results/guild_results.tsv y results/diamond_results.tsv
```

### Ejemplo 3: Solo GUILD con red personalizada
```bash
python scripts/network_propagation.py \
  -n data/network_guild.txt \
  -g 2023 5230 3099 \
  -a guild \
  -o resultados_guild
```

## Parámetros de los Algoritmos

### GUILD (modificable en el código)
```python
def guild(G, genes_semilla, r=0.5, iteraciones=100, epsilon=1e-6):
```

- **r**: Probabilidad de restart (0-1)
  - Valor bajo (0.2): Explora más lejos
  - Valor alto (0.8): Se queda cerca de semillas
  - **Default: 0.5** (balance)

- **iteraciones**: Máximo de iteraciones (default: 100)
- **epsilon**: Criterio de convergencia (default: 1e-6)

### DIAMOnD (modificable en el código)
```python
def diamond(G, genes_semilla, top_k=200):
```

- **top_k**: Número de genes candidatos a identificar
  - **Default: 200**
  - Ajustar según tamaño de la red

## Funcionalidades Avanzadas

### Agregar Genes Faltantes Automáticamente

Si tus genes semilla no están en la red, el script intenta:

1. Convertir símbolos → Entrez IDs
2. Buscar interacciones en STRING
3. Agregar automáticamente a la red

```bash
python scripts/network_propagation.py -n data/network_diamond.txt -g ENO1 PGK1 HK2
```

```
Genes no encontrados: ENO1, PGK1, HK2
Intentando conversión automática a Entrez IDs...
   ENO1 → 2023
   PGK1 → 5230
   HK2 → 3099
   Buscando interacciones en STRING...
   Encontradas 6 interacciones
   Agregadas 3 interacciones a data/network_diamond.txt
```

## Solución de Problemas

### Error: "Ningún gen semilla en la red"

**Causa**: Los genes no están en la red o hay incompatibilidad de IDs

**Solución**:
```bash
# Verificar qué tipo de IDs usa tu red
head -5 data/network_guild.txt

# Si usa Entrez IDs, convertir genes:
python scripts/network_propagation.py -g ENO1 PGK1 HK2  # Se convierten automáticamente

# Si usa símbolos, usar red STRING:
python scripts/network_propagation.py -n data/string_network_filtered_hugo-400.tsv -g ENO1 PGK1 HK2
```

### Error: "Formato no reconocido"

**Causa**: El archivo de red tiene formato no estándar

**Solución**: Verificar que siga uno de los formatos soportados:
- GUILD: `nodo1 peso nodo2`
- DIAMOnD: `nodo1,nodo2`
- STRING: TSV con encabezados

### DIAMOnD devuelve p-valores altos (>0.05)

**Causas comunes**:
1. Muy pocos genes semilla (< 5)
2. Red muy densa (muchas conexiones)
3. Genes semilla son "hubs" (muy conectados)

**Soluciones**:
```bash
# 1. Usar más genes semilla (10-20 recomendado)
python scripts/network_propagation.py -g ENO1 PGK1 HK2 GAPDH TPI1 ALDOA PKM LDHA

# 2. Confiar en GUILD para estos casos
python scripts/network_propagation.py -g ENO1 PGK1 HK2 -a guild
```

### Error: "mygene no instalado"

**Solución**:
```bash
pip install mygene
```
## Casos de Uso

### 1. Identificar genes candidatos para enfermedad
```bash
# Genes conocidos de una enfermedad
python scripts/network_propagation.py -g BRCA1 BRCA2 TP53 ATM -a both
```

### 2. Expandir vía metabólica
```bash
# Genes de glucólisis → identificar genes relacionados
python scripts/network_propagation.py -g ENO1 PGK1 HK2 GAPDH TPI1 -a guild
```

### 3. Análisis de red específica de tejido
```bash
# Usar red personalizada de tejido específico
python scripts/network_propagation.py -n data/brain_network.txt -s data/genes.txt
```

## Referencias

- **STRING Database**: https://string-db.org/
- **GUILD Algorithm**: https://github.com/emreg00/guild
- **DIAMOnD Algorithm**: https://github.com/dinaghiassian/DIAMOnD
- **MyGene.info**: https://mygene.info/

## Licencia

Este proyecto está bajo la Licencia MIT.

## Autor

Hugo Salas Calderón - [hugosalascalderon@gmail.com](mailto:hugosalascalderon@gmail.com)

---

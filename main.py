from clustering import hierarchical_clustering, compute_cophenetic_coefficient, draw_dendrograms
import networkx as nx
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import substitution_matrices
from scipy.cluster.hierarchy import dendrogram
from Bio.SubsMat import MatrixInfo as matlist
import io
import matplotlib.pyplot as plt
import numpy as np
from flask import Flask, render_template, request, jsonify, url_for
import base64
import matplotlib
from BLOSUM62 import BLOSUM62
matplotlib.use('Agg')





app = Flask(__name__)


#####################################################
############# Analisis de secuencias ################
#####################################################

def es_adn(secuencia):
    return all(nuc in 'ATCG' for nuc in secuencia)


def es_arn(secuencia):
    return all(nuc in 'AUCG' for nuc in secuencia)


def es_proteina(secuencia):
    aminoacidos = 'ACDEFGHIKLMNPQRSTVWY'
    return all(aa in aminoacidos for aa in secuencia)


def transcribir_adn_a_arn(secuencia):
    return secuencia.replace('T', 'U')


@app.route('/secuencias')
def secuencias():
    return render_template('secuencias.html')


@app.route('/procesar', methods=['POST'])
def procesar():
    secuencia = request.form['secuencia'].upper()
    tipo = ''
    transcripcion = ''
    imagen_url = ''

    if es_adn(secuencia):
        tipo = 'ADN'
        transcripcion = transcribir_adn_a_arn(secuencia)
        imagen_url = url_for('static', filename='images/adn.png')
    elif es_arn(secuencia):
        tipo = 'ARN'
        imagen_url = url_for('static', filename='images/arn.png')
    elif es_proteina(secuencia):
        tipo = 'Proteína'
        imagen_url = url_for('static', filename='images/proteina.png')
    else:
        tipo = 'Desconocido'
        imagen_url = url_for('static', filename='images/desconocido.png')

    resultado = {
        'tipo': tipo,
        'cantidad': len(secuencia),
        'transcripcion': transcripcion if tipo == 'ADN' else 'N/A',
        'imagen_url': imagen_url
    }

    return jsonify(resultado)


#####################################################
######### Matriz puntos ############
#####################################################

def create_dot_matrix(seq1, seq2):
    matrix = np.zeros((len(seq1), len(seq2)))
    for i, char1 in enumerate(seq1):
        for j, char2 in enumerate(seq2):
            if char1 == char2:
                matrix[i, j] = 1
    return matrix


def plot_dot_matrix(matrix, seq1, seq2):
    fig, ax = plt.subplots(figsize=(8, 6))
    rows, cols = np.where(matrix == 1)
    ax.scatter(cols, rows, c='blue', marker='o',
               edgecolor='black', s=100, alpha=0.7)
    ax.set_xticks(np.arange(len(seq2)))
    ax.set_yticks(np.arange(len(seq1)))
    ax.set_xticklabels(list(seq2), fontsize=12)
    ax.set_yticklabels(list(seq1), fontsize=12)
    ax.set_xlabel('Secuencia 2', fontsize=14)
    ax.set_ylabel('Secuencia 1', fontsize=14)
    ax.set_title('Matriz de Puntos', fontsize=16)
    ax.grid(True, linestyle='--', color='gray', alpha=0.7)
    ax.invert_yaxis()
    plt.xticks(rotation=45, ha='right')
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plt.close(fig)
    return buf


@app.route('/matrizPuntos')
def matrizPuntos():
    return render_template('matrizPuntos.html')


@app.route('/dotmatrix', methods=['POST'])
def dotmatrix():
    seq1 = request.form['seq1']
    seq2 = request.form['seq2']
    matrix = create_dot_matrix(seq1, seq2)
    buf = plot_dot_matrix(matrix, seq1, seq2)
    image = base64.b64encode(buf.getvalue()).decode('ascii')
    return jsonify({'image': image})


#####################################################
######### Neddleman-Wunch ############
#####################################################


@app.route('/nw')
def nw():
    return render_template('neddleman-wunch/index.html')


#####################################################
######### Alineamiento proteínas (usar Blosum o Pam) - FASTA ############
#####################################################

blosum62 = substitution_matrices.load("BLOSUM62")


def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5)
    best_alignment = alignments[0]
    return best_alignment


@app.route('/blosum', methods=['GET', 'POST'])
def blosum():
    if request.method == 'POST':
        seq1 = request.form['sequence1']
        seq2 = request.form['sequence2']
        alignment = align_sequences(seq1, seq2)
        return render_template('blosum.html', alignment=alignment, seq1=seq1, seq2=seq2)
    return render_template('blosum.html')

#####################################################
######### Alineamiento proteínas (usar Blosum o Pam) - FASTA ############
#####################################################

#####################################################
######### START Alignment ############
#####################################################


# Penalizaciones
gap = -2
match = 1
mismatch = -1

# Extraer la matriz de alineación y máximo score de dos alineaciones (Needleman Wunsch parcial)


def max_score(seq_1, seq_2):
    left, up = (seq_1, seq_2) if len(seq_1) > len(seq_2) else (seq_2, seq_1)

    matrix = [[0 for j in range(len(up)+1)] for i in range(len(left)+1)]
    for i in range(len(matrix)):
        matrix[i][0] = i*gap
    for j in range(len(matrix[0])):
        matrix[0][j] = j*gap

    for i in range(1, len(matrix)):
        for j in range(1, len(matrix[0])):
            matrix[i][j] = max(matrix[i-1][j-1]+(match if left[i-1] == up[j-1]
                               else mismatch), matrix[i-1][j]+gap, matrix[i][j-1]+gap)
    return matrix, matrix[-1][-1]

# Extraer la posición de la cadena estrella


def get_star_pos(sequences):
    pair_scores = [[0 for j in range(len(sequences)+1)]
                   for i in range(len(sequences)+1)]
    pair_matrices = [[[] for j in range(len(sequences)+1)]
                     for i in range(len(sequences)+1)]
    best_score = float("-inf")
    best_score_pos = -1

    for i in range(len(sequences)):
        for j in range(i, len(sequences)):
            if i == j:
                continue
            matrix, score = max_score(sequences[i], sequences[j])
            pair_matrices[i][j] = matrix
            pair_matrices[j][i] = matrix

            pair_scores[i][j] = score
            pair_scores[j][i] = score
        pair_scores[i][len(sequences)] = sum(pair_scores[i])
        pair_scores[len(sequences)][i] = pair_scores[i][len(sequences)]

        if pair_scores[i][len(sequences)] > best_score:
            best_score = pair_scores[i][len(sequences)]
            best_score_pos = i

    scores = [pair_scores[i][len(sequences)] for i in range(len(sequences))]

    return best_score_pos, best_score, pair_scores, pair_matrices, scores

# Obtener una alineación óptima de 2 secuencias en base a su matriz de alineación


def align_sequences(seq_1, seq_2, matrix):
    AlignmentA = ""
    AlignmentB = ""
    i = len(seq_1)
    j = len(seq_2)

    while i > 0 and j > 0:
        Score = matrix[i][j]
        ScoreDiag = matrix[i-1][j-1] if i > 0 and j > 0 else float('-inf')
        ScoreUp = matrix[i][j-1] if j > 0 else float('-inf')
        ScoreLeft = matrix[i-1][j] if i > 0 else float('-inf')

        if Score == ScoreDiag + (match if seq_1[i-1] == seq_2[j-1] else mismatch):
            AlignmentA = seq_1[i-1] + AlignmentA
            AlignmentB = seq_2[j-1] + AlignmentB
            i -= 1
            j -= 1
        elif Score == ScoreLeft + gap:
            AlignmentA = seq_1[i-1] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i -= 1
        else:  # Score == ScoreUp + gap
            AlignmentA = "-" + AlignmentA
            AlignmentB = seq_2[j-1] + AlignmentB
            j -= 1

    while i > 0:
        AlignmentA = seq_1[i-1] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i -= 1

    while j > 0:
        AlignmentA = "-" + AlignmentA
        AlignmentB = seq_2[j-1] + AlignmentB
        j -= 1

    return AlignmentA, AlignmentB


def insert_substrig(pos, seq, subseq='-'):
    return seq[:pos] + subseq + seq[pos:]

# Ejecutar alineación estrella y obtener alineación múltiple


def star_alignment(sequences):
    star_pos, best_score, pair_scores, pair_matrices, scores = get_star_pos(
        sequences)

    alignments = ['' for i in range(len(sequences))]
    star_alignments = []

    for i in range(len(sequences)):
        if star_pos == i:
            continue
        star_ali, i_ali = align_sequences(
            sequences[star_pos], sequences[i], pair_matrices[star_pos][i])
        star_alignments.append(
            (star_ali, i_ali, i, pair_matrices[star_pos][i][-1][-1]))
        if star_ali != alignments[star_pos] or len(star_ali) != len(alignments[star_pos]):
            for k in range(len(alignments[star_pos])):
                if alignments[star_pos][k] == '-':
                    prev = 0
                    while prev < i-1:
                        if len(alignments[prev]) < len(alignments[star_pos]):
                            alignments[prev] = insert_substrig(
                                k, alignments[prev])
                        prev += 1
                    i_ali = insert_substrig(k, i_ali)
                    star_ali = insert_substrig(k, star_ali)
        alignments[star_pos] = star_ali
        alignments[i] = i_ali

    return star_pos, best_score, pair_scores, star_alignments, alignments, scores


@app.route('/estrella', methods=['GET', 'POST'])
def estrella():
    alignments = None
    star_alignments = None
    pair_scores = None
    sequences = None
    best_score = None
    star_pos = None
    scores = None
    if request.method == 'POST':
        sequences = request.form.get('sequences').strip().split('\n')
        sequences = [seq.strip() for seq in sequences]
        star_pos, best_score, pair_scores, star_alignments, alignments, scores = star_alignment(
            sequences)
    return render_template('estrella.html', sequences=sequences, star_pos=star_pos, best_score=best_score, pair_scores=pair_scores, star_alignments=star_alignments, alignments=alignments, scores=scores)

#####################################################
######### BLOSUM ############
#####################################################

#####################################################
######### ARBOL ENRAIZADO ############
#####################################################


#####################################################
######### CLUSTERIZACIÓN ############
#####################################################
@app.route('/clusterizacion')
def clusterizacion():
    return render_template('form.html')


@app.route('/results', methods=['POST'])
def results():
    data = request.form['matrix']
    method = request.form['method']

    # Convert the input string into a matrix
    distance_matrix = np.array(
        [[float(num) for num in line.split()] for line in data.strip().split('\n')])

    if method == 'Compare':
        # Compare all methods
        methods = ['Minimum', 'Maximum', 'Average']
        results = {}
        clusters_list = []
        distances_list = []

        for method in methods:
            clustering_result, matrices, min_distances, cophenetic_matrix = hierarchical_clustering(
                distance_matrix, method)
            ccc = compute_cophenetic_coefficient(
                distance_matrix, cophenetic_matrix)
            results[method] = (clustering_result, matrices,
                               min_distances, cophenetic_matrix, ccc)
            clusters_list.append(clustering_result)
            distances_list.append(min_distances)

        best_method = max(results, key=lambda m: results[m][4])
        best_ccc = results[best_method][4]

        # Draw the dendrogram
        draw_dendrograms(clusters_list, distances_list, methods,
                         save_filename='static/images/combined_dendrograms.png')

        return render_template('compare_results.html',
                               results=results,
                               best_method=best_method,
                               best_ccc=best_ccc)
    else:
        # Perform clustering for the selected method
        clustering_result, matrices, min_distances, cophenetic_matrix = hierarchical_clustering(
            distance_matrix, method)
        ccc = compute_cophenetic_coefficient(
            distance_matrix, cophenetic_matrix)

        # Draw the dendrogram
        draw_dendrograms([clustering_result], [min_distances], [
                         method], save_filename='static/images/combined_dendrograms.png')

        return render_template('results.html',
                               clusters=clustering_result,
                               min_distances=min_distances,
                               cophenetic_matrix=cophenetic_matrix,
                               ccc=ccc,
                               method=method)


#####################################################
######### SECONDARY STRUCTURE ############
#####################################################

def compute_energy_matrix(sequence, alpha):
    length = len(sequence)
    energy_matrix = [[0] * length for _ in range(length)]

    for gap in range(1, length):
        for i in range(length - gap):
            j = i + gap
            min_energy = float('inf')
            min_energy = min(min_energy, energy_matrix[i + 1][j])
            min_energy = min(min_energy, energy_matrix[i][j - 1])
            if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA', 'GU', 'UG']:
                min_energy = min(
                    min_energy, energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]])
            for k in range(i + 1, j):
                min_energy = min(
                    min_energy, energy_matrix[i][k - 1] + energy_matrix[k][j])

            energy_matrix[i][j] = min_energy

    return energy_matrix


def predict_secondary_structure(sequence, alpha):
    length = len(sequence)
    energy_matrix = compute_energy_matrix(sequence, alpha)
    score = energy_matrix[0][length - 1]

    return energy_matrix, score


def traceback(sequence, energy_matrix, alpha):
    length = len(sequence)
    traceback_pairs = ["."] * length
    paired_positions = []
    structures = []

    def traceback_helper(i, j):
        nonlocal traceback_pairs, paired_positions, structures

        if i >= j:
            return

        if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA', 'GU', 'UG'] and energy_matrix[i][j] == energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]]:
            traceback_pairs[i] = "("
            traceback_pairs[j] = ")"
            paired_positions.append((i + 1, j + 1))
            traceback_helper(i + 1, j - 1)
            return

        if energy_matrix[i][j] == energy_matrix[i + 1][j]:
            traceback_helper(i + 1, j)
            return

        if energy_matrix[i][j] == energy_matrix[i][j - 1]:
            traceback_helper(i, j - 1)
            return

        for k in range(i + 1, j):
            if energy_matrix[i][j] == energy_matrix[i][k - 1] + energy_matrix[k][j]:
                traceback_helper(i, k - 1)
                traceback_helper(k, j)
                return

    traceback_helper(0, length - 1)

    # Identificar estructuras: bucles internos, tallos, protuberancias
    visited = [False] * length
    for i, pair in enumerate(traceback_pairs):
        if pair == "(" and not visited[i]:
            j = traceback_pairs.index(")", i + 1)
            visited[i] = visited[j] = True
            if j - i > 1:
                if all(traceback_pairs[x] == "." for x in range(i + 1, j)):
                    structures.append(f"Lazo: Bucle interno: {i + 1}-{j + 1}")
                else:
                    structures.append(f"Tallo: {i + 1}-{j + 1}")
            else:
                if i == 0 or j == length - 1:
                    structures.append(
                        f"Lazo: Base no emparejada: {i + 1}-{j + 1}")
                else:
                    structures.append(f"Bulbo: {i + 1}-{j + 1}")

    return traceback_pairs, paired_positions, structures


def plot_structure(sequence, pairs):
    plt.figure(figsize=(max(4, len(sequence) / 8), max(4, len(sequence) / 8)))

    G = nx.Graph()
    node_colors = []
    base_color_map = {'A': '#FF0000', 'U': '#00FF00', 'C': '#0000FF', 'G': '#FFFF00'}

    for idx, base in enumerate(sequence):
        G.add_node(idx + 1, base=base, color=base_color_map[base])
        node_colors.append(base_color_map[base])

    for idx in range(1, len(sequence)):
        G.add_edge(idx, idx + 1, color="#808080", style="-")

    for pair in pairs:
        G.add_edge(pair[0], pair[1], color="#FF5733", style="--")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, edgecolors='black', linewidths=1.5)
    nx.draw_networkx_labels(G, pos, labels={i: f"{i} ({G.nodes[i]['base']})" for i in G.nodes()}, font_size=10, font_color="black")

    edges = G.edges()
    edge_colors = [G[u][v]['color'] for u, v in edges]
    edge_styles = [G[u][v]['style'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, style=edge_styles, width=2)

    # Agregar la leyenda fuera de la imagen
    plt.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=base)
        for base, color in base_color_map.items()
    ] + [plt.Line2D([0], [0], color="#FF5733", lw=2, linestyle='--', label='Par emparejado')],
    loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)

    plt.axis('off')

    # Guardar la imagen en un objeto de bytes en memoria
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    plt.close()

    # Convertir la imagen en base64
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')

    return image_base64



@app.route('/secondary', methods=['GET', 'POST'])
def secondary():
    if request.method == 'POST':
        sequence = request.form['sequence'].upper()
        alpha_dict = {"CG": -1, "GC": -1, "AU": -
                      1, "UA": -1, "GU": -1, "UG": -1}

        energy_matrix, min_energy_score = predict_secondary_structure(
            sequence, alpha_dict)
        traceback_structure, paired_positions, structures = traceback(
            sequence, energy_matrix, alpha_dict)

        image_base64 = plot_structure(sequence, paired_positions)

        return render_template('secondary.html', sequence=sequence, energy_matrix=energy_matrix, min_energy_score=min_energy_score, paired_positions=paired_positions, structures=structures, image_base64=image_base64)

    return render_template('secondary.html')


#####################################################
######### BLOSUM ############
#####################################################

def algoritmo_alineamiento_local(sequence_1, sequence_2, gap_penalty):
    len_1 = len(sequence_1)
    len_2 = len(sequence_2)

    array = np.zeros(shape=(len_1 + 1, len_2 + 1))

    max_score = 0
    max_i = 0
    max_j = 0

    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            array[i, j] = max(0,
                              array[i - 1, j - 1] +
                              BLOSUM62[sequence_1[i-1]][sequence_2[j-1]],
                              array[i - 1, j] + gap_penalty,
                              array[i, j - 1] + gap_penalty)
            if array[i, j] > max_score:
                max_score = array[i, j]
                max_i = i
                max_j = j

    return array, max_score, max_i, max_j

def algoritmo_alineamiento_global(sequence_1, sequence_2, gap_penalty):
    len_1 = len(sequence_1)
    len_2 = len(sequence_2)

    array = np.zeros(shape=(len_1 + 1, len_2 + 1))

    for i in range(1, len_1 + 1):
        array[i, 0] = gap_penalty * i
    for j in range(1, len_2 + 1):
        array[0, j] = gap_penalty * j

    for i in range(1, len_1 + 1):
        for j in range(1, len_2 + 1):
            array[i, j] = max(array[i - 1, j - 1] + BLOSUM62[sequence_1[i-1]][sequence_2[j-1]],
                              array[i - 1, j] + gap_penalty,
                              array[i, j - 1] + gap_penalty)

    max_score = array[len_1, len_2]
    return array, max_score, len_1, len_2

def traceback_global(array, sequence_1, sequence_2, gap_penalty):
    aligned_seq_1 = ""
    aligned_seq_2 = ""
    i = len(sequence_1)
    j = len(sequence_2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and array[i, j] == array[i - 1, j - 1] + BLOSUM62[sequence_1[i-1]][sequence_2[j-1]]:
            aligned_seq_1 = sequence_1[i - 1] + aligned_seq_1
            aligned_seq_2 = sequence_2[j - 1] + aligned_seq_2
            i -= 1
            j -= 1
        elif i > 0 and array[i, j] == array[i - 1, j] + gap_penalty:
            aligned_seq_1 = sequence_1[i - 1] + aligned_seq_1
            aligned_seq_2 = '-' + aligned_seq_2
            i -= 1
        else:
            aligned_seq_1 = '-' + aligned_seq_1
            aligned_seq_2 = sequence_2[j - 1] + aligned_seq_2
            j -= 1

    return aligned_seq_1, aligned_seq_2

def traceback_local(array, sequence_1, sequence_2, max_i, max_j, gap_penalty):
    aligned_seq_1 = ""
    aligned_seq_2 = ""
    i = max_i
    j = max_j

    while array[i, j] != 0:
        if array[i, j] == array[i - 1, j - 1] + BLOSUM62[sequence_1[i-1]][sequence_2[j-1]]:
            aligned_seq_1 = sequence_1[i - 1] + aligned_seq_1
            aligned_seq_2 = sequence_2[j - 1] + aligned_seq_2
            i -= 1
            j -= 1
        elif array[i, j] == array[i - 1, j] + gap_penalty:
            aligned_seq_1 = sequence_1[i - 1] + aligned_seq_1
            aligned_seq_2 = '-' + aligned_seq_2
            i -= 1
        elif array[i, j] == array[i, j - 1] + gap_penalty:
            aligned_seq_1 = '-' + aligned_seq_1
            aligned_seq_2 = sequence_2[j - 1] + aligned_seq_2
            j -= 1

    start_i = i + 1  # Convert to 1-based index
    start_j = j + 1  # Convert to 1-based index

    return aligned_seq_1, aligned_seq_2, start_i, start_j

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/align', methods=['POST'])
def align():
    sequence_1 = request.form['sequence1']
    sequence_2 = request.form['sequence2']
    gap_penalty = int(request.form['gap_penalty'])
    alignment_type = request.form['alignment_type']

    if alignment_type == 'local':
        array, max_score, max_i, max_j = algoritmo_alineamiento_local(sequence_1, sequence_2, gap_penalty)
        aligned_seq_1, aligned_seq_2, start_i, start_j = traceback_local(array, sequence_1, sequence_2, max_i, max_j, gap_penalty)
        return render_template('index.html', alignment_type=alignment_type, sequence_1=sequence_1, sequence_2=sequence_2,
                               aligned_seq_1=aligned_seq_1, aligned_seq_2=aligned_seq_2, max_score=max_score, start_i=start_i, start_j=start_j,
                               score_matrix=array.tolist())
    else:
        array, max_score, _, _ = algoritmo_alineamiento_global(sequence_1, sequence_2, gap_penalty)
        aligned_seq_1, aligned_seq_2 = traceback_global(array, sequence_1, sequence_2, gap_penalty)
        return render_template('index.html', alignment_type=alignment_type, sequence_1=sequence_1, sequence_2=sequence_2,
                               aligned_seq_1=aligned_seq_1, aligned_seq_2=aligned_seq_2, max_score=max_score,
                               score_matrix=array.tolist())

#####################################################
######### Final ############
#####################################################


if __name__ == '__main__':
    app.run(debug=True)

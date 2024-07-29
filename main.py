import base64
from flask import Flask, render_template, request, jsonify, url_for

import numpy as np
import matplotlib.pyplot as plt
import io
import base64
from Bio.SubsMat import MatrixInfo as matlist


# test
from Bio.Align import substitution_matrices
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

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
######### Final ############
#####################################################


if __name__ == '__main__':
    app.run(debug=True)

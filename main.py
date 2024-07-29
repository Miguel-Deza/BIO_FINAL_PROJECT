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

# BLOSUM62
BLOSUM62 = {
    'A': {'A':  4, 'R': -1, 'N': -2, 'D': -2, 'C':  0, 'Q': -1, 'E': -1, 'G':  0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  0, 'W': -3, 'Y': -2, 'V':  0},
    'R': {'A': -1, 'R':  5, 'N':  0, 'D': -2, 'C': -3, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R':  0, 'N':  6, 'D':  1, 'C': -3, 'Q':  0, 'E':  0, 'G':  0, 'H':  1, 'I': -3, 'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P': -2, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N':  1, 'D':  6, 'C': -3, 'Q':  0, 'E':  2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C':  9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R':  1, 'N':  0, 'D':  0, 'C': -3, 'Q':  5, 'E':  2, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  1, 'M':  0, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R':  0, 'N':  0, 'D':  2, 'C': -4, 'Q':  2, 'E':  5, 'G': -2, 'H':  0, 'I': -3, 'L': -3, 'K':  1, 'M': -2, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A':  0, 'R': -2, 'N':  0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G':  6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S':  0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R':  0, 'N':  1, 'D': -1, 'C': -3, 'Q':  0, 'E':  0, 'G': -2, 'H':  8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y':  2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I':  4, 'L':  2, 'K': -3, 'M':  1, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V':  3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I':  2, 'L':  4, 'K': -2, 'M':  2, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V':  1},
    'K': {'A': -1, 'R':  2, 'N':  0, 'D': -1, 'C': -3, 'Q':  1, 'E':  1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K':  5, 'M': -1, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q':  0, 'E': -2, 'G': -3, 'H': -2, 'I':  1, 'L':  2, 'K': -1, 'M':  5, 'F':  0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V':  1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I':  0, 'L':  0, 'K': -3, 'M':  0, 'F':  6, 'P': -4, 'S': -2, 'T': -2, 'W':  1, 'Y':  3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P':  7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A':  1, 'R': -1, 'N':  1, 'D':  0, 'C': -1, 'Q':  0, 'E':  0, 'G':  0, 'H': -1, 'I': -2, 'L': -2, 'K':  0, 'M': -1, 'F': -2, 'P': -1, 'S':  4, 'T':  1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A':  0, 'R': -1, 'N':  0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  5, 'W': -2, 'Y': -2, 'V':  0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F':  1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y':  2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H':  2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F':  3, 'P': -3, 'S': -2, 'T': -2, 'W':  2, 'Y':  7, 'V': -1},
    'V': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I':  3, 'L':  1, 'K': -2, 'M':  1, 'F': -1, 'P': -2, 'S': -2, 'T':  0, 'W': -3, 'Y': -1, 'V':  4}
}

def algoritmo_alineamiento_local_blosum(sequence_1, sequence_2, gap_penalty):
    len_1 = len(sequence_1)
    len_2 = len(sequence_2)

    array = np.zeros(shape=(len_1 + 1, len_2 + 1))

    max_score = 0
    max_i = 0
    max_j = 0

    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            array[i, j] = max(0,
                              array[i - 1, j - 1] + BLOSUM62[sequence_1[i-1]][sequence_2[j-1]],
                              array[i - 1, j] + gap_penalty,
                              array[i, j - 1] + gap_penalty)
            if array[i, j] > max_score:
                max_score = array[i, j]
                max_i = i
                max_j = j

    return array, max_score, max_i, max_j

def traceback_local_blosum(array, sequence_1, sequence_2, max_i, max_j, gap_penalty):
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


# ENDPOINTS 

@app.route('/')
def index():
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
        score = align_sequences(seq1, seq2)
        return render_template('blosum.html', score=score, seq1=seq1, seq2=seq2)
    return render_template('blosum.html')

#####################################################
######### Alineamiento proteínas (usar Blosum o Pam) - FASTA ############
#####################################################


#####################################################
######### Final ############
#####################################################


if __name__ == '__main__':
    app.run(debug=True)

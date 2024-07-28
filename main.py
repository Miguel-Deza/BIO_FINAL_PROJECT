import base64
from flask import Flask, render_template, request, jsonify, url_for

import numpy as np
import matplotlib.pyplot as plt
import io
import base64
from Bio.SubsMat import MatrixInfo as matlist

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
######### Final ############
#####################################################


if __name__ == '__main__':
    app.run(debug=True)

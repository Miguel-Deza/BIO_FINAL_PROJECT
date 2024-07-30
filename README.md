# 🧬 PROYECTO FINAL BIOINFORMÁTICA VISUALIZADOR 🧬

¡Bienvenido al proyecto final de Bioinformática! 🎉 Este proyecto es un visualizador que te ayudará a entender y analizar datos bioinformáticos de manera interactiva y visual. 🌟

## 📁 Estructura del Proyecto

```plaintext
BIO_FINAL_PROJECT
├─ .vscode
│  └─ settings.json
├─ BLOSUM62.py
├─ clustering.py
├─ main.py
├─ README.md
├─ requirements.txt
├─ static
│  ├─ GridBuilder.js
│  ├─ GridBuilder_w.js
│  ├─ images
│  │  ├─ clusterizacion_combined_dendrograms.png
│  │  ├─ combined_dendrograms.png
│  │  ├─ d.png
│  │  ├─ ds.png
│  │  ├─ dsu.png
│  │  ├─ du.png
│  │  ├─ s.png
│  │  ├─ su.png
│  │  └─ u.png
│  ├─ structure.png
│  ├─ style.css
│  └─ style_w.css
├─ templates
│  ├─ blosum.html
│  ├─ cantidadelementos.html
│  ├─ compare_results.html
│  ├─ enraizado.html
│  ├─ estrella.html
│  ├─ form.html
│  ├─ home.html
│  ├─ input.html
│  ├─ is_possible_transcription.html
│  ├─ matrizPuntos.html
│  ├─ neddleman-wunch
│  │  └─ index.html
│  ├─ nj.html
│  ├─ results.html
│  ├─ secondary.html
│  ├─ secuencias.html
│  ├─ smith-waterman
│  │  ├─ images
│  │  │  ├─ d.png
│  │  │  ├─ ds.png
│  │  │  ├─ dsu.png
│  │  │  ├─ du.png
│  │  │  ├─ s.png
│  │  │  ├─ su.png
│  │  │  └─ u.png
│  │  └─ index.html
│  ├─ upgma.html
│  └─ whatis.html
├─ test.py
└─ __pycache__
   ├─ BLOSUM62.cpython-312.pyc
   └─ clustering.cpython-312.pyc

```
## 🚀 Instalación

Para instalar y ejecutar este proyecto, sigue estos pasos:

1. **Clona el repositorio**:
    ```bash
    git clone https://github.com/Miguel-Deza/BIO_FINAL_PROJECT.git
    ```
2. **Navega al directorio del proyecto**:
    ```bash
    cd BIO_FINAL_PROJECT
    ```
3. **Crea un entorno virtual** (recomendado):
    ```bash
    python -m venv venv
    ```
4. **Activa el entorno virtual**:
    - En Windows:
        ```bash
        venv\Scripts\activate
        ```
    - En macOS y Linux:
        ```bash
        source venv/bin/activate
        ```
5. **Instala las dependencias**:
    ```bash
    pip install -r requirements.txt
    ```

## 🧩 Uso

Para ejecutar el proyecto, usa el siguiente comando:

```bash
python main.py
```

## 📚 Funcionalidades

- **Matrices de BLOSUM62**: Genera y visualiza matrices de sustitución de aminoácidos. 🔬
- **Algoritmos de Clustering**: Realiza clustering jerárquico para análisis filogenéticos. 🌳
- **Algoritmos de Alineamiento**: Implementación de algoritmos Needleman-Wunsch y Smith-Waterman para alineamiento de secuencias. 🔄
- **Visualización Interactiva**: Interfaz web intuitiva para visualizar los resultados de los análisis. 🖥️

## 🤝 Contribuciones

¡Las contribuciones son bienvenidas! Si deseas contribuir, por favor sigue estos pasos:

1. Haz un fork del proyecto.
2. Crea una nueva rama (`git checkout -b feature-nueva-funcionalidad`).
3. Realiza tus cambios y haz commit (`git commit -m 'Agrega nueva funcionalidad'`).
4. Haz push a la rama (`git push origin feature-nueva-funcionalidad`).
5. Abre un Pull Request.
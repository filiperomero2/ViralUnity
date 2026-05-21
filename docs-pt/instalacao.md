# Instalação

O ViralUnity é um pacote Python que executa workflows Snakemake. As dependências da ferramenta estão documentadas no arquivo `environment.yml` e as dependências específicas de cada workflow estão em `viralunity/scripts/envs/`.

## Clone e crie o ambiente

```bash
git clone https://github.com/InstitutoTodosPelaSaude/ViralUnity.git
cd ViralUnity/
conda env create -f environment.yml
conda activate viralunity
```

Ou com **micromamba** (recomendado no macOS com Apple Silicon):

```bash
micromamba env create -f environment.yml --platform osx-64
micromamba activate viralunity
```

```{warning}
No macOS com Apple Silicon (M1 ou posterior), o ambiente `viralunity/scripts/envs/clair3.yaml` pode não ser instalado corretamente devido a restrições de compatibilidade nas dependências do clair3.
```

## Verifique a instalação

```bash
viralunity --version
viralunity --help
```

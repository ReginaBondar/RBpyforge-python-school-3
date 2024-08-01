FROM continuumio/miniconda3
RUN conda install conda-forge::rdkit
WORKDIR /code
COPY ./main.py /code
COPY ./requirements.txt /app
RUN pip install -r /code/requirements.txt
EXPOSE 8000
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]

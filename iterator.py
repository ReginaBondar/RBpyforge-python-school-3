class IteratorMolecules:
    def __init__(self, sequence, molecule, stop=1):
        self._sequence = sequence
        self._stop = stop
        self._mol = molecule
        self._index = 0  # счетчик для перебора элементов последовательности
        self._counter = 0  # количество возвращенных значений


    def __iter__(self):
        return self

    def __next__(self):
        if self._counter < self._stop and self._index < len(self._sequence):
            if self._sequence[self._index].HasSubstructMatch(self._mol):
                item = self._index
                self._index += 1
                self._counter += 1
                return item
            else:
                self._index += 1
        else:
            raise StopIteration

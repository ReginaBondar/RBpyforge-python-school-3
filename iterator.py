class IteratorMolecules:
    def __init__(self, sequence, stop=1):
        self._sequence = sequence
        self._stop = stop
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._index < self._stop:
            molecule = self._sequence[self._index]
            self._index += 1
            return molecule
        else:
            raise StopIteration

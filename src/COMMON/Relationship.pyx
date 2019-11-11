cdef class Relationship:
    cdef float _in, _out
    cdef int num_edges_to, num_edges_from

    def __init__(self, float _in, int num_edges_to, float _out, int num_edges_from):
        self._in = _in
        self.num_edges_to = num_edges_to
        self._out = _out
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self._in, self.num_edges_to, self._out, self.num_edges_from)





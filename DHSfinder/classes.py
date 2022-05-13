#Defines classes for sequence features

class SeqUnit:

    def __init__(self, contig, sequence, position, sequence_type, family, identifier):

        self.sequence = sequence
        self.position = int(position)
        self.contig = contig
        self.sequence_type = sequence_type
        self.family = family
        self.identifier = identifier

    # if the sequence type is nuc then the position is the actual position in the genome of the kmer and if the sequence type is protein then the distance
    # is the number of genes

    def can_it_join(self, other_sequnit, min_space, max_distance):

        return (self.contig == other_sequnit.contig and abs(self.position - other_sequnit.position) < max_distance and abs(self.position - other_sequnit.position) > min_space)


class AcessoryRegion:

    def __init__(self, first_bound_family, second_bound_family):

        self.first_bound_family = first_bound_family
        self.second_bound_family = second_bound_family
        self.regions = []
        self.score = 0


class RegionID:

    def __init__(self, first_bound, second_bound):

        self.first_bound = first_bound
        self.second_bound = second_bound

    def __eq__(self, other):

        return (self.first_bound == other.first_bound or self.first_bound == other.second_bound) and (self.second_bound == other.first_bound or self.second_bound == other.second_bound)

    def __hash__(self):

        return hash(self.first_bound) ^ hash(self.second_bound)


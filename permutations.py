from gridclasses import VectorSet


__all__ = ['Permutation', 'PegPerm']

start_index = 1  # for display purposes only

# todo make all this stuff just tuples, maybe zip together entries/signs
# todo check the combination stuff (vector bases) and enumeration stuff


class Permutation:
    """ Object representing permutations. Permutations are automatically
    standardized to 0 through n-1 on creation, which is used extensively
    later on. """

    def __init__(self, entries, signs=None, b=None):  # TODO is b the right symbol

        if isinstance(entries, Permutation):  # clone constructor
            self.entries = entries.entries
            self.signs = entries.signs
            self.b = entries.b
            return

        assert len(set(entries)) == len(entries), 'make sure elements are distinct!'
        assert len(signs) == len(entries), 'length of signs and entries should match!'

        if isinstance(entries, tuple) or isinstance(entries, list):
            sorted_entries = sorted(entries)
            self.entries = tuple(map(lambda e: sorted_entries.index(e), entries))
        else:
            # TODO raise error
            pass

        if signs:
            self.signs = tuple(signs)
        else:
            # TODO raise error
            pass

        if b:
            self.b = b
        else:
            # TODO raise error
            pass

    def __repr__(self):
        return ''.join([str(a+start_index) + '\''*(b+1) + '\t' for a, b in zip(self.entries, self.signs)])

    def __eq__(self, other):
        if isinstance(other, Permutation):
            return self.entries == other.entries and self.signs == other.signs
        else:
            return False

    def __hash__(self):
        return hash(self.entries) + 2*hash(self.signs)

    def __len__(self):
        return len(self.entries)

    def delete(self, idx) -> 'Permutation':
        return Permutation(self.entries[:idx] + self.entries[idx+1:], self.signs[:idx] + self.signs[idx+1:], self.b)

    def plot(self):  # TODO update
        """ Draws a plot of the given Permutations. """
        n = len(self)
        array = [[' ' for i in range(n)] for j in range(n)]
        for i in range(n):
            array[self[i]][i] = '*'
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        # return s
        print(s)

    def fills(self):  # TODO check this
        pp = PegPerm(self, [0]*len(self))
        pp, vectorset = pp.clean()  # TODO this should call new method "reduce" that does what their "compact" did
        # the * unpacks the set, so we zip together vectors
        # the absolute value is magic: if there is a dot by itself, the compact
        # function puts a 0 in the vector, then the -1 and abs turns it to a 1
        vector = [abs(max(t)-1) for t in zip(*vectorset.basis)]
        # pp expanded by the vector should match the original permutation
        assert pp.inflate(vector) == self
        return pp, vector


class PegPerm:
    """ A permutation with 'pegged' entries, representing classes of permutations.
    Pegs can be +, ., or -. Object represented by a sign tuple and a permutation. """

    def __init__(self, perm, decorators=None):

        if isinstance(perm, PegPerm):  # clone constructor
            self.perm = Permutation(perm.perm)
            self.decorators = perm.decorators
            return

        if isinstance(perm, Permutation) and isinstance(decorators, list):  # TODO decorators could be tuple
            self.perm = Permutation(perm)
        else:
            raise Exception('Bad inputs for PegPerm constructor')

        self.perm = perm
        self.decorators = tuple(decorators)  # TODO MAKE SURE THEY ARE +/- 1

    def __len__(self):
        return len(self.decorators)

    def __repr__(self):
        signsymbols = ['(.)', '(+)', '(-)']
        return ''.join([str(a+start_index) + '\''*(b+1) + signsymbols[c] + '\t'
                        for a, b, c in zip(self.perm.entries, self.perm.signs, self.decorators)])

    def __hash__(self):
        return hash(self.perm) + 2*hash(self.decorators)

    def __eq__(self, other):
        if isinstance(other, PegPerm):
            return self.decorators == other.decorators and self.perm == other.perm
        else:
            return False

    def plot(self):  # TODO FINISH
        n = self.__len__()
        array = [[' ' for i in range(n)] for j in range(n)]
        signlist = ['.', '+', '-']
        for i in range(n):
            array[self.perm[i]][i] = str(signlist[self.signs[i]])
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        # return s
        print(s)

    def minimum_fill(self):
        newentries = []
        newsigns = []

        for entry, sign, decorator in zip(self.perm.entries, self.perm.signs, self.decorators):
            newentries.append(3*entry)
            newsigns.append(sign)
            if decorator != 0:
                newentries.append(3*entry + decorator)
                newsigns.append(sign)

        return Permutation(newentries, newsigns)

    def inflate(self, vector) -> Permutation:
        assert len(vector) == len(self)

        new_entries = []
        new_signs = []

        multiplier = 2*max(vector)

        for entry, sign, decorator, number in zip(self.perm.entries, self.perm.signs, self.decorators, vector):
            for i in range(number):
                new_entries.append(multiplier*entry + decorator*i)
                new_signs.append(sign)
        return Permutation(new_entries, new_signs, self.perm.b)

    def delete(self, m):
        return PegPerm(self.perm.delete(m), self.decorators[:m]+self.decorators[m+1:])

    def dot(self, m):
        return PegPerm(self.perm, self.decorators[:m]+(0,)+self.decorators[m+1:])

    def dotall(self):
        p = self.perm
        n = len(self)
        S = set()

        for i in range(1, 2**n):

            # builds binary vector from i
            k = i
            v = [0 for x in range(n)]
            for m in range(n):
                if k >= 2**(n-m-1):
                    v[m]=1
                    k = k-2**(n-m-1)

            # dots permutation using the vector
            q = PegPerm(self)
            for j in range(n):
                if v[j] == 1:
                    q = q.dot(j)  # TODO this could certainly be faster if we construct all these perms 1 dot at a time

            # adds q to the set
            S.add(q)

        return S

    def complete(self) -> set['PegPerm']:  # TODO could this be faster?
        n = len(self)
        L = [set() for i in range(n + 1)]
        L[n] = self.dotall()
        L[n].add(self)
        while n >= 1:
            for q in L[n]:
                for i in range(len(q)):
                    L[n-1].add(q.delete(i))
            n -= 1
        return set().union(*L)

    def is_compact(self) -> bool:  # TODO write cleaner
        for i in range(len(self)-1):
            if self.perm.signs[i] == self.perm.signs[i+1]:
                a = (self.perm.entries[i]+2) * self.decorators[i]
                b = (self.perm.entries[i+1]+2) * self.decorators[i+1]

                if b-a == 1 or b == 0:
                    return False
        return True

    def clean(self) -> tuple['PegPerm', VectorSet]:  # intended to work on PegPerms that are known to be compact

        entries = []
        signs = []
        decorators = []
        vector = []  # essentially the "maximum" vector we could inflate by, where 0 represents infinity

        idx = 0
        while idx < len(self):
            if self.decorators[idx] == 0:
                val = self.signs  # TODO

            else:
                entries.append(self.perm.entries[idx])
                signs.append(self.perm.signs[idx])
                decorators.append(self.decorators[idx])
                vector.append(0)
                idx += 1

            val = entries[idx + 1] - entries[idx]
            # if a bond, checks if both signs are dots. if not, continues
            if abs(val) == 1:
                runlength = 1
                while (decorators[idx + runlength - 1] == 0 and
                       decorators[idx + runlength] == 0 and
                       entries[idx + runlength] - entries[idx + runlength - 1] == val and
                       ):
                    runlength += 1
                    # stop if we hit the end of the permutation
                    if idx + runlength >= len(signs):
                        break
                # if we have a run, delete those entries and add to vector
                if runlength > 1:
                    for i in range(runlength - 1):
                        perm = perm.delete(idx + 1)
                        del signs[idx + 1]
                    del vector[idx + 1 : idx + runlength]
                    vector[idx] = runlength + 1
                    signs[idx] = val
            idx += 1

        # create avoidance basis
        basis = set()
        for idx, val in enumerate(vector):
            if val:
                new_vector = [0 for i in range(len(entries))]
                new_vector[idx] = val
                basis.add(tuple(new_vector))

        # minimum vector is just the minimum vector that fills our new perm (idk if this is natural but that's what the paper describes)
        minimum = [1 + abs(i) for i in signs]

        return PegPerm(Permutation(entries, signs, self.perm.b), decorators), VectorSet(minimum, basis)

    def split_blocks(self, *splits):  # helper
        splits = list(splits)
        splits.sort()
        splits.reverse()
        n = len(splits)

        new_decorators = list(self.decorators)
        mult = 3  # TODO this could prob be 1
        new_entries = [mult*n*x for x in self.perm.entries]
        new_signs = list(self.perm.signs)

        spacer = n
        for split in splits:
            sign, val = self.decorators[split], self.perm.entries[split]
            new_decorators.insert(split+1, sign)
            new_entries.insert(split+1, mult*n*val + spacer*sign)
            new_signs.insert(split+1, self.perm.signs[split])
            spacer -= 1
        return new_entries, new_signs, new_decorators

    def block_transpose(self):  # TODO later
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    signs, perm = self.split_blocks(i, j, k)
                    signs = signs[:i+1] + signs[j+2:k+3] + signs[i+1:j+2] + signs[k+3:]
                    perm = perm[:i+1] + perm[j+2:k+3] + perm[i+1:j+2] + perm[k+3:]
                    S.add(PegPerm(signs, perm))
        return S

    def cut_paste(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one cut paste
                move from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    signs, perm = self.split_blocks(i, j, k)
                    # break up sign into pieces
                    s1 = signs[:i+1]; s2 = signs[j+2:k+3]
                    s3 = signs[i+1:j+2]; s4 = signs[k+3:]
                    # break up perm into pieces
                    p1 = perm[:i+1]; p2 = perm[j+2:k+3]
                    p3 = perm[i+1:j+2]; p4 = perm[k+3:]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
                    # reverse piece two, add it to set
                    s2.reverse(); p2.reverse()
                    s2 = [-sign for sign in s2]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
                    # undo the reversal, then do the same for piece 3
                    s2.reverse(); p2.reverse()
                    s2 = [-sign for sign in s2]
                    s3.reverse(); p3.reverse()
                    s3 = [-sign for sign in s3]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
        return S

    def block_interchange(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one block
                interchange from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    for m in range(k, n):
                        if self.signs[m] == 0:
                            continue
                        s, p = self.split_blocks(i, j, k, m)
                        s = s[:i+1] + s[k+3:m+4] + s[j+2:k+3] + s[i+1:j+2] + s[m+4:]
                        p = p[:i+1] + p[k+3:m+4] + p[j+2:k+3] + p[i+1:j+2] + p[m+4:]
                        S.add(PegPerm(s, p))
        return S

    def block_reversal(self):
        """ Returns the set of PegPerms which can result from one block
        reversal from the current one. """
        S = set([])
        n = len(self)

        for i in range(n):
            if self.decorators[i] == 0:
                continue
            for j in range(i, n):
                if self.decorators[j] == 0:
                    continue

                perm, signs, decorators = self.split_blocks(i, j)
                b = self.perm.b

                p_start, p_mid, p_end = perm[:i+1], perm[i+1:j+2], perm[j+2:]
                s_start, s_mid, s_end = signs[:i+1], [(s+1) % b for s in signs[i+1:j+2]], signs[j+2:]
                d_start, d_mid, d_end = decorators[:i+1], [-d for d in decorators[i+1:j+2]], decorators[j+2:]

                p_mid.reverse()
                s_mid.reverse()
                d_mid.reverse()

                S.add(PegPerm(Permutation(p_start + p_end, s_start + s_end, b), d_start + d_end))
        return S

    def prefix_transpose(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one block
                reversal from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                signs, perm = self.split_blocks(i,j)
                s_start, s_mid, s_end = signs[:i+1], signs[i+1:j+2], signs[j+2:]
                p_start, p_mid, p_end = perm[:i+1], perm[i+1:j+2], perm[j+2:]
                signs = s_mid + s_start + s_end
                perm = p_mid + p_start + p_end
                S.add(PegPerm(signs, perm))
        return S

    def prefix_reversal(self):
        """ Returns the set of PegPerms which can result from one prefix
        reversal from the current one. """
        S = set()
        n = len(self)

        for i in range(n):
            if self.decorators[i] == 0:
                continue

            perm, signs, decorators = self.split_blocks(i)
            b = self.perm.b

            p_start, p_end = perm[:i+1], perm[i+1:]
            s_start, s_end = [(s+1) % b for s in signs[:i+1]], signs[i+1:]
            d_start, d_end = [-d for d in decorators[:i+1]], decorators[i+1:]

            p_start.reverse()
            s_start.reverse()
            d_start.reverse()

            S.add(PegPerm(Permutation(p_start+p_end, s_start+s_end, b), d_start+d_end))
        return S

    @staticmethod
    def sort_pr(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.prefix_reversal())
            S = S_
        return S

    @staticmethod
    def sort_pt(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.prefix_transpose())
            S = S_
        return S

    @staticmethod
    def sort_br(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_reversal())
            S = S_
        return S

    @staticmethod
    def sort_bt(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_transpose())
            S = S_
        return S

    @staticmethod
    def sort_bi(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_interchange())
            S = S_
        return S

    @staticmethod
    def sort_cp(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.cut_paste())
            S = S_
        return S

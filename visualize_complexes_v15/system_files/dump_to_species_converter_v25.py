from struct import unpack

# code parameters
agents, agent, left, right, comp = ['A', 'T', 'C', 'G'], 2, 5, 7, 9
nuc = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}


def dump_to_species(inputfile, des_path, filename):
    def uni(dumpfile):
        with open(dumpfile, 'rb') as f:
            data = f.read()
            size = int(len(data)/8)
            f.close()

        values = [int(unpack('<d', data[i * 8:i * 8 + 8])[0]) for i in range(0, size)]
        lenn = values.index(values[0] + 1)
        li = (values[i:i+lenn] for i in range(0, len(values), lenn))
        return li

    r = sorted(list(uni(inputfile)), key=lambda s: int(s[1]))

    def transform(a):
        le = '5!' + str(a[1]) if a[1] else '5'
        ri = '3!' + str(a[2]) if a[2] else '3'
        w = 'W!' + str(a[3]) if a[3] else 'W'

        return 'N(b~{},{},{},{})'.format(a[0], le, ri, w)

    global ssdna, ssdna_1, vq, n, c
    ssdna, ssdna_1, vcx, vq, n, c = [], [], [], {}, 1, 0
    comp_n = [i[1] for i in r]

    for e in r:
        if e[1] not in vcx:
            vcx.append(e[1])
            ssdna.append(ssdna_1) if ssdna_1 != [] else None
            ssdna_1, vq, n, c = [], {}, 1, comp_n.count(e[1]) + 1

        def set_state(cl):
            global n, c
            sets = []
            le, ri, w = cl[left], cl[right], cl[comp]

            if le != -1:
                if le in vq:
                    item = vq[le][1]
                    sets.append(item)

                elif le not in vq:
                    sets.append(n)
                    n += 1

            elif le == -1:
                sets.append(None)

            if ri != -1:
                if ri in vq:
                    item = vq[ri][0]
                    sets.append(item)

                elif ri not in vq:
                    sets.append(n)
                    n += 1

            elif ri == -1:
                sets.append(None)

            if w != -1:
                if w in vq:
                    item = vq[w][2]
                    sets.append(item)

                elif w not in vq:
                    sets.append(c)
                    c += 1

            elif w == -1:
                sets.append(None)

            return sets

        items = set_state(e)
        bngl_syn = transform([nuc[e[2]], items[0], items[1], items[2]])
        ssdna_1.append(bngl_syn)
        vq.update({e[0]: [items[0], items[1], items[2]]})
        ssdna.append(ssdna_1) if e == r[-1] else None

    new, nn = [], []
    for c in ssdna:
        if c not in new:
            new.append(c)
            nn.append(1)

        elif c in new:
            ind = new.index(c)
            nn[ind] += 1

    with open(des_path, 'w') as f:
        f.write("%s" % '# Species list generated by VDNA for file: ' + "'" + filename + "'" + '\n')
        f.write("%s" % '\n')
        for c, cc in zip(new, nn):
            item = '.'.join(c) + '  ' + str(cc) + '\n'
            f.write("%s" % item)
    f.close()

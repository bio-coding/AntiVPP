import sys
from PyQt5 import uic, QtWidgets
from sklearn.externals import joblib

qtCreatorFile = "antivpp.ui"  # Name of the file here

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def partialDictRoundedSum(d, keys):
    return round(sum([d[k] for k in keys]), 3)


def getResultsSeq(sequence, rfc):
    paste_seq = str(sequence)
    seq_size = len(paste_seq+str(0.000001))

    kyte_doolittle = {'A': 1.80, 'C': 2.50, 'D': -3.50, 'E': -3.50, 'F': 2.80,
                      'G': -0.40, 'H': -3.20, 'I': 4.50, 'K': -3.90, 'L': 3.80,
                      'M': 1.90, 'N': -3.50, 'P': -1.60, 'Q': -3.50, 'R': -4.50,
                      'S': -0.80, 'T': -0.70, 'V': 4.20, 'W': -0.90, 'Y': -1.30}

    molecular_weigth = {'A': 89.09, 'C': 121.15, 'D': 133.10, 'E': 147.13, 'F': 165.19,
                        'G': 75.07, 'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17,
                        'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
                        'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.24, 'Y': 181.19}

    net_charge = {'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0,
                  'G': 0, 'H': 0, 'I': 0, 'K': 1, 'L': 0,
                  'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1,
                  'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

    net_hydrogen = {'A': 0, 'C': 0, 'D': 1, 'E': 1, 'F': 0,
                    'G': 0, 'H': 1, 'I': 0, 'K': 2, 'L': 0,
                    'M': 0, 'N': 2, 'P': 0, 'Q': 2, 'R': 4,
                    'S': 1, 'T': 1, 'V': 0, 'W': 1, 'Y': 1}

    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    aa_counts = {k: paste_seq.count(k) for k in aa_list}

    aa_perc = {k: round(aa_counts[k]/seq_size, 3) for k in aa_list}

    # PROPERTIES Q-P

    aliphatic = partialDictRoundedSum(aa_perc, ['I', 'V', 'L'])

    negative_charged = partialDictRoundedSum(aa_perc, ['D', 'E'])

    total_charged = partialDictRoundedSum(aa_perc, ['D', 'E', 'K', 'H', 'R'])

    aromatic = partialDictRoundedSum(aa_perc, ['F', 'H', 'W', 'Y'])

    polar = partialDictRoundedSum(aa_perc, ['D', 'E', 'R', 'K', 'Q', 'N'])

    neutral = partialDictRoundedSum(aa_perc,
                                    ['A', 'G', 'H', 'P', 'S', 'T', 'Y'])

    hydrophobic = partialDictRoundedSum(aa_perc,
                                        ['C', 'F', 'I', 'L', 'M', 'V', 'W'])

    positive_charged = partialDictRoundedSum(aa_perc, ['K', 'R', 'H'])

    tiny = partialDictRoundedSum(aa_perc, ['A', 'C', 'D', 'G', 'S', 'T'])

    small = partialDictRoundedSum(aa_perc,
                                  ['E', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'V'])

    large = partialDictRoundedSum(aa_perc, ['F', 'R', 'W', 'Y'])

    # SCALES

    kyleD = round(
        sum(
            [aa_counts[k]*kyte_doolittle[k] for k in aa_list]
        )/seq_size, 3)

    molW = round(sum([aa_counts[k]*molecular_weigth[k] for k in aa_list]), 3)

    netCharge = sum([aa_counts[k]*net_charge[k] for k in aa_list])

    netH = round(sum([aa_counts[k]*net_hydrogen[k] for k in aa_list]), 3)

    result = rfc.predict([
        [netH, netCharge, molW, kyleD] +
        [v for v in aa_perc.values()] +
        [tiny, small, large, aliphatic, aromatic,
         total_charged, negative_charged, positive_charged,
         polar, neutral, hydrophobic]])

    return str(result)


def getResultsFile(filename):
    rfc = joblib.load('modelo_entrenado_2.pkl')

    with open(filename) as fp:
        print('Name\tPredicted_Antiviral\tSequence\n')
        for name, seq in read_fasta(fp):
            result = getResultsSeq(seq, rfc)
            print('{}\t{}\t{}'.format(
                name.replace('>', ''),
                result.replace('[', '').replace(']', '').strip(),
                seq
            ))


class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.setFixedSize(self.size())  # Dimensiones fijas
        # Esto es para ordenar que cuando se presione vaya a calculo
        self.bottonpred.clicked.connect(self.calculation)

    def calculation(self):
        from sklearn.externals import joblib
        rfc = joblib.load('modelo_entrenado_2.pkl')

        paste_seq = str(self.pasteseq.toPlainText())

        seq_size = len(paste_seq+str(0.000001))

        kyte_doolittle = {'A': 1.80, 'C': 2.50, 'D': -3.50, 'E': -3.50, 'F': 2.80,
                          'G': -0.40, 'H': -3.20, 'I': 4.50, 'K': -3.90, 'L': 3.80,
                          'M': 1.90, 'N': -3.50, 'P': -1.60, 'Q': -3.50, 'R': -4.50,
                          'S': -0.80, 'T': -0.70, 'V': 4.20, 'W': -0.90, 'Y': -1.30}

        molecular_weigth = {'A': 89.09, 'C': 121.15, 'D': 133.10, 'E': 147.13, 'F': 165.19,
                            'G': 75.07, 'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17,
                            'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
                            'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.24, 'Y': 181.19}

        net_charge = {'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0,
                      'G': 0, 'H': 0, 'I': 0, 'K': 1, 'L': 0,
                      'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1,
                      'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

        net_hydrogen = {'A': 0, 'C': 0, 'D': 1, 'E': 1, 'F': 0,
                        'G': 0, 'H': 1, 'I': 0, 'K': 2, 'L': 0,
                        'M': 0, 'N': 2, 'P': 0, 'Q': 2, 'R': 4,
                        'S': 1, 'T': 1, 'V': 0, 'W': 1, 'Y': 1}

        aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        aa_counts = {k: paste_seq.count(k) for k in aa_list}

        aa_perc = {k: round(aa_counts[k]/seq_size, 3) for k in aa_list}

        # PROPERTIES Q-P

        aliphatic = partialDictRoundedSum(aa_perc, ['I', 'V', 'L'])

        negative_charged = partialDictRoundedSum(aa_perc, ['D', 'E'])

        total_charged = partialDictRoundedSum(
            aa_perc, ['D', 'E', 'K', 'H', 'R'])

        aromatic = partialDictRoundedSum(aa_perc, ['F', 'H', 'W', 'Y'])

        polar = partialDictRoundedSum(aa_perc, ['D', 'E', 'R', 'K', 'Q', 'N'])

        neutral = partialDictRoundedSum(aa_perc,
                                        ['A', 'G', 'H', 'P', 'S', 'T', 'Y'])

        hydrophobic = partialDictRoundedSum(aa_perc,
                                            ['C', 'F', 'I', 'L', 'M', 'V', 'W'])

        positive_charged = partialDictRoundedSum(aa_perc, ['K', 'R', 'H'])

        tiny = partialDictRoundedSum(aa_perc, ['A', 'C', 'D', 'G', 'S', 'T'])

        small = partialDictRoundedSum(aa_perc,
                                      ['E', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'V'])

        large = partialDictRoundedSum(aa_perc, ['F', 'R', 'W', 'Y'])

        # SCALES

        kyleD = round(
            sum(
                [aa_counts[k]*kyte_doolittle[k] for k in aa_list]
            )/seq_size, 3)

        molW = round(sum([aa_counts[k]*molecular_weigth[k]
                          for k in aa_list]), 3)

        netCharge = sum([aa_counts[k]*net_charge[k] for k in aa_list])

        netH = round(sum([aa_counts[k]*net_hydrogen[k] for k in aa_list]), 3)

        result = rfc.predict([
            [netH, netCharge, molW, kyleD] +
            [v for v in aa_perc.values()] +
            [tiny, small, large, aliphatic, aromatic,
             total_charged, negative_charged, positive_charged,
             polar, neutral, hydrophobic]])

        result = "Probable: " + str(result)

        self.textpred.setText(str(result))
        self.textpred1.setText(str(aliphatic))
        self.textpred2.setText(str(negative_charged))
        self.textpred3.setText(str(aromatic))
        self.textpred4.setText(str(polar))
        self.textpred5.setText(str(neutral))
        self.textpred6.setText(str(hydrophobic))
        self.textpred7.setText(str(positive_charged))
        self.textpred8.setText(str(tiny))
        self.textpred9.setText(str(small))
        self.textpred10.setText(str(large))
        self.textpred11.setText(str(kyleD))
        self.textpred12.setText(str(molW))
        self.textpred13.setText(str(netCharge))
        self.textpred14.setText(str(netH))
        self.textpred15.setText(str(total_charged))
        self.textrelat.setText(' , '.join(
            [k+': '+str(aa_perc[k]) for k in aa_list]))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        app = QtWidgets.QApplication(sys.argv)
        window = MyApp()
        window.show()
        sys.exit(app.exec_())
    else:
        if sys.argv[-1] == '--help':
            print(
                "Usage:\n"
                "For UI version:\n"
                "\tpython antivpp.py\n"
                "For getting results from fasta file:\n"
                "\tpython antivpp.py [fasta filename] > [tsv filename]\n"
                "\nUse relative paths\n")
        else:
            print(getResultsFile(sys.argv[-1]))

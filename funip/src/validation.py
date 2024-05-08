from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from io import StringIO

# from .logger import Mes

# validate if string is good sequence
def validate_seq(seqstring):

    if seqstring.startswith(">"):
        tmp = StringIO(seqstring)

        try:
            seqlist = list(SeqIO.parse(tmp, fasta))
            return "seqrecord", seqlist
        except:
            logging.warning(f"Invalid seqrecord {seqstring}")
            return "invalid", None

    else:
        seqstring = seqstring.replace(" ", "").replace("\n", "")

        try:
            seq = Seq(seqstring)
            seqlist = [
                SeqRecord(
                    seq,
                    id="input",
                    description="tmp",
                    annotations={"molecule_type": "DNA"},
                )
            ]
            return "seqrecord", seqlist
        except:
            logging.warning(f"Invalid sequence {seqstring}")
            return "invalid", None

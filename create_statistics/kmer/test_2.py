from minipy import KMer, minimize
from Bio import SeqIO
import torch
from tokenizers import SentencePieceBPETokenizer
from transformers import PreTrainedTokenizerFast


EG_SEQ ="AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC" \
          "AACCTTGGACTACGATCGGGGGRACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC"

KMER_LEN = 6 #10 #29
WIN_LEN = 5 #9

#minimizers = minimize(EG_SEQ, KMER_LEN, WIN_LEN)
#for m in minimizers:
#    print(f'{m.value()}, {m.position()}, {m.strand()}')

"""
Traverse a chromosome for minimizers.
Store all minimizers in a dict with the count of appearence
--> hope for enough RAM
choose the 1Mio minimizers with highest occurance rate
check in the chromosome and on a 2nd chromosome, how much percent of the minimizers is covered by the dict
do this for different k-mer sizes
"""
def tokenize_sequence(seq):
    tokens = minimize(seq, KMER_LEN, WIN_LEN)
    token_list = []
    for t in tokens:
        token_list.append(t.value())
    return token_list

"""sequences = []
for record in SeqIO.parse("../data/mio_snippets_chr17/overfit/s10/reads_chr17_10_seqreq.fasta", "fasta"):
    sequences.append(str(record.seq))
#for record in SeqIO.parse("chr21.fasta", "fasta"):
#    test_sequence = str(record.seq)
lens = []
for s in sequences:
    token_sequence = tokenize_sequence(s)
    lens.append(len(token_sequence))
print(f'max length: {max(lens)} ,avg length: {sum(lens)/len(lens)}')"""


####################################Sentence Piece####################################################################

"""sequences = []
for record in SeqIO.parse("../data/chm13.v1.1.fasta", "fasta"):
    sequences.append(str(record.seq))
    break
text = sequences[0]

special_tokens = ["<s>", "<pad>", "</s>", "<unk>", "<cls>", "<sep>", "<mask>"]
tk_tokenizer = SentencePieceBPETokenizer()
print("init")
tk_tokenizer.train_from_iterator(
    text,
    vocab_size=4096,
    min_frequency=2,
    show_progress=True,
    special_tokens=special_tokens
)
tk_tokenizer.save('tokens_sentencepiece')
# convert
tokenizer = PreTrainedTokenizerFast(tokenizer_object=tk_tokenizer, special_tokens=special_tokens)
tokenizer.bos_token = "<s>"
tokenizer.bos_token_id = tk_tokenizer.token_to_id("<s>")
tokenizer.pad_token = "<pad>"
tokenizer.pad_token_id = tk_tokenizer.token_to_id("<pad>")
tokenizer.eos_token = "</s>"
tokenizer.eos_token_id = tk_tokenizer.token_to_id("</s>")
tokenizer.unk_token = "<unk>"
tokenizer.unk_token_id = tk_tokenizer.token_to_id("<unk>")
tokenizer.cls_token = "<cls>"
tokenizer.cls_token_id = tk_tokenizer.token_to_id("<cls>")
tokenizer.sep_token = "<sep>"
tokenizer.sep_token_id = tk_tokenizer.token_to_id("<sep>")
tokenizer.mask_token = "<mask>"
tokenizer.mask_token_id = tk_tokenizer.token_to_id("<mask>")
# and save for later!
tokenizer.save_pretrained("./path/to/transformers/version/")
####################################Sentence Piece####################################################################


print("finished")
#print(f'{m.value()}, {m.position()}, {m.strand()}')"""
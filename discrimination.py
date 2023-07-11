#!/usr/bin/python
###############################################################################
# Program expects as input...
#
# --fasta   Fasta file that has identifies formated such:
#
#           ><unique-id> <genus> <species>
#           <..sequence..>
#           ><unique-id> <genus> <species>
#           <..sequence..>
#
#           The full fasta file is used to build the the BLAST database. Only
#           sequences for species that have >= 2 sequences will be used within
#           the descrimination analysis.
#
# --meta    Excel meta-data that has the following column header names
#           (at least). Note that "order" is optional but the rest are
#           mandatory.
#
#           Accepted_Name, [order], family, genus, species
#
#           The Accepted_Name is used to match against the <species> <genus>
#           from the fasta file.
###############################################################################

import argparse
import json
import logging

from pandas import DataFrame, read_excel, read_csv
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from subprocess import check_output
from collections import defaultdict


# Setup logging to output to file and to console...
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%m-%d %H:%M',
    filename='discrimination-analysis-debug.log',
    filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
logging.getLogger('').addHandler(console)
logger = logging.getLogger(__name__)


class FastaSpeciesException(Exception):
    pass


class MetaFileException(Exception):
    def __init__(self, fname, columns):
        self.message = f'\nMissing column(s) in {fname}\n'
        self.message += '\n'.join(columns)
        super().__init__(self.message)


class CheckFile(argparse.Action):
    '''
    Action for argparse to check that the file exists.
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        if not Path(values).is_file():
            parser.error(f"The file {values} does not exist!")
        else:
            setattr(namespace, self.dest, Path(values))


def read_meta(meta_fname: Path) -> DataFrame:
    '''
    Checks the meta file has all the requested headers and if so returns the
    the data frame.

    If not an exception is raised with an indication of the misssing headers.
    '''
    meta_taxon_map = {
        'Accepted_Name': 'a',
        'family': 'f',
        'genus': 'g',
        'species': 's'
    }
    meta_key_map = {v: k for k, v in meta_taxon_map.items()}

    meta_data_df = read_excel(meta_fname)
    taxon = ''

    for k, v in meta_taxon_map.items():
        taxon += (v if k in meta_data_df else '')

    if 'afgs' not in taxon:
        # Missin headers!
        missing = [v for k, v in meta_key_map.items() if k not in taxon]
        raise MetaFileException(meta_fname, missing)

    return meta_data_df


def read_fasta(fasta_fname: str) -> list[SeqRecord]:
    '''
    Reads the FASTA file and returns a list of sequences
    '''
    with open(fasta_fname, 'r') as fd:
        return [seq for seq in SeqIO.parse(fd, "fasta")]


def annotate_with_taxonomy(
        sequences: list[SeqRecord],
        metadata: DataFrame) -> list[SeqRecord]:
    '''
    We run through the sequences making sure that the species exists
    within the metadata and adding in the taxonomy for our distrimination.
    We return the anonated sequences or raise an exception.
    '''
    annotated = []

    for sequence in sequences:
        # The description will be <id> <genus> <species>
        # We want to make sure that the <genus> <species>
        # is listed in the "Accepted_Name" of the metadata
        parts = sequence.description.split()
        speices = ' '.join(parts[1:])

        # Grab the taxonomy information
        result = metadata.loc[metadata['Accepted_Name'] == speices]

        if len(result) != 1:
            raise FastaSpeciesException(f'"{speices}" not in Metadata')

        # We found it so annotate the description to use later in
        # the BLAST result. Note that we hyphante the species, if not
        # already to make things easier when matching in the BLAST.
        annotation = {
            'Accepted_Name': result['Accepted_Name'].to_string(index=False),
            'family': result['family'].to_string(index=False),
            'genus': result['genus'].to_string(index=False),
            'species': '-'.join(
                result['species'].to_string(index=False).split())
        }

        # Add in order if it exists...
        if 'order' in metadata:
            annotation['order'] = result['order'].to_string(index=False)

        # We quote the annotation so it is treated as a single string
        # Note that "hopefully" a forward slash will not be found in a
        # taxonomy name.
        sequence.description = (
            f'/{json.dumps(annotation, separators=(",", ":"))}/')
        annotated.append(sequence)

    return annotated


def save_sequences(sequences: list[SeqRecord]) -> Path:
    '''
    Saves the sequence to file and returns the path to that file, the sequence
    is saved unaligned.
    '''
    fname = Path('blast_discrimination_fasta.fas')

    with open(fname, 'w') as fd:
        for seq in sequences:
            fd.write(f'>{seq.id} {seq.description}\n')
            fd.write(f'{seq.seq.replace("-", "")}\n')

    return fname


def make_blast_database(seq_fname: Path) -> str:
    '''
    We are assuming we have access to the "makeblastdb" and it is on
    the path. We make a directory and save off the sequences
    and make the BLAST directory and return the directory name.

    Note on Ubuntu: sudo apt-get install ncbi-blast+
    '''
    makeblastdb_bin = Path('makeblastdb')
    databasename = 'blastdb'
    cmd = [
        f'{makeblastdb_bin}',
        "-in",
        f'{seq_fname}',
        "-out",
        f'{databasename}',
        "-dbtype",
        "nucl",
        "-title",
        "temp_blastdb",
        "-parse_seqids"
    ]
    check_output(cmd)

    return databasename


def blastn(fname_query, dbname='blastdb', out='output.csv') -> list[dict]:
    '''
    Run a query on the blastdb name given using the filename for the
    query sequence.

    Return a list of dict of IDs and stitles that match the top bitscores.
    Note that the stitles will have the JSON taxonomy information.
    '''
    cmd = [
        'blastn',
        '-db',
        dbname,
        '-negative_seqidlist',
        fname_query,
        '-query',
        fname_query,
        '-out',
        out,
        '-outfmt',
        '10 std score stitle',
        '-max_target_seqs',
        '100'
    ]
    check_output(cmd)

    headers = [
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore',
        'score',
        'stitle',
    ]

    df = read_csv(out, header=None, names=headers, quotechar='/')
    logger.debug('<====== Top 20 results =====>')
    logger.debug(f'{df.head(20)}')

    # Grab the top bit scores
    logger.debug('<====== Top bitscore results =====>')
    results = df[df['bitscore'] == df['bitscore'].max()]
    logger.debug(f'{results}')

    top_matches = []
    for _, row in results.iterrows():
        top_matches.append(
            {
                'id': row['sseqid'],
                'taxonomy': json.loads(row['stitle'])
            }
        )

    return top_matches


def blast_sequence(
        sequence: SeqRecord,
        dbname: Path) -> list[dict]:
    '''
    Given a sequence we blast it and pass back the results. To run the
    blast we have to save the sequence to file first and then BLAST it
    against the dabasebase.
    '''
    blast_sequence_fname = 'query.fas'

    logger.debug('<===== BLAST SEQUENCE ======>')
    logger.debug('Blasting: %s' % (sequence.id))
    logger.debug('Blasting: %s' % (sequence.description))

    with open(blast_sequence_fname, 'w') as fd:
        fd.write(f'>{sequence.id}\n')
        fd.write(f'{sequence.seq.replace("-", "")}\n')

    return blastn(blast_sequence_fname)


def discriminate(sequence: SeqRecord, results: list[dict]) -> list:
    '''
    Returns a list of the discrimination results for this sequence based upon
    the blast results. This is sutable for loading stright into a DataFrame
    '''
    if len(results) == 0:
        logger.warning(
            f'!!!! No Blast Results for: "{sequence.id}" !!!!')

    taxonomy = json.loads(sequence.description.strip('/'))

    discrim = {
        'sid': sequence.id,
        'taxa': taxonomy['Accepted_Name'],
        'order': 'N/A',
        'family': taxonomy['family'],
        'genus': taxonomy['genus'],
        'species': taxonomy['species'],
        'order_uniq': True,
        'family_uniq': True,
        'genus_uniq': True,
        'species_uniq': True
    }

    for result in results:
        if "order" in result['taxonomy']:
            if taxonomy['order'] != result['taxonomy']['order']:
                discrim['order_uniq'] = False
                discrim['family_uniq'] = False
                discrim['genus_uniq'] = False
                discrim['species_uniq'] = False
        if taxonomy['family'] != result['taxonomy']['family']:
            discrim['family_uniq'] = False
            discrim['genus_uniq'] = False
            discrim['species_uniq'] = False
        if taxonomy['genus'] != result['taxonomy']['genus']:
            discrim['genus_uniq'] = False
            discrim['species_uniq'] = False
        if taxonomy['species'] != result['taxonomy']['species']:
            discrim['species_uniq'] = False

    return [
            discrim['sid'],
            discrim['taxa'],
            discrim['order'],
            discrim['family'],
            discrim['genus'],
            discrim['species'],
            discrim['order_uniq'],
            discrim['family_uniq'],
            discrim['genus_uniq'],
            discrim['species_uniq']
   ]


def analysis(sequences: list[SeqRecord], dbname: str) -> list[list]:
    '''
    Using the database that has been created and the sequences we
    blast each sequence that has a >= 2 for the species and calculate
    the discrimination against the taxonomy.

    Returns a list of lists of discrimination results.
    '''

    discrimination = []

    # Collect by accepted name
    sequences_by_species = defaultdict(list)
    for sequence in sequences:
        taxonomy = json.loads(sequence.description.strip('/'))
        sequences_by_species[taxonomy['Accepted_Name']].append(sequence)

    # Where we have more than one sequence for an accepted name we blast it.
    for species, seq_list in sequences_by_species.items():
        if len(seq_list) > 1:
            logger.debug('<====== Processing Species ======>')
            logger.debug(f'--> {species}')
            for sequence in seq_list:
                logger.debug('<====== Processing Sequence ======>')
                logger.debug(f'--> {sequence.id}')

                discrimination.append(
                    discriminate(
                        sequence,
                        blast_sequence(sequence, dbname)))
        else:
            logger.debug('<====== Not Processing ======>')
            logger.debug(f'"{species}" as only 1 sequence')

    return discrimination


def save_discrimination(
        discrim: list[dict],
        fname='discrimination_result.xlsx'):
    '''
    Saves the discrimination results as an Excel spreadsheet.
    '''
    df = DataFrame(
        discrim,
        columns=['ID',
                 'Taxa',
                 'Order',
                 'Family',
                 'Genus',
                 'Species',
                 'Order Match',
                 'Family Match',
                 'Genus Match',
                 'Species Match'])

    df.to_excel(fname, engine='xlsxwriter')
    logger.info(f'Results are in: {fname}')


def run_discrimination(meta_fname: Path, fasta_fname: Path):
    '''
    To run the discrimination we must first build the BLAST database and
    pick out the data we need from the meta file. Once done we can then
    BLAST each sequence in turn against the database and output the result
    as a new excel file showing the descimination for a given sequence
    to the taxonomy defined i.e., family, genus, species.
    '''
    annotated_sequences = annotate_with_taxonomy(
        read_fasta(fasta_fname), read_meta(meta_fname))

    save_discrimination(
        analysis(
            annotated_sequences,
            make_blast_database(save_sequences(annotated_sequences))))


if __name__ == '__main__':
    argp = argparse.ArgumentParser(
        description='Run discrimination analytics on FASTA through BLAST')
    argp.add_argument(
        '-f',
        '--fasta',
        action=CheckFile,
        help='Path to Fasta file',
        required=True)
    argp.add_argument(
        '-m',
        '--meta',
        action=CheckFile,
        help='Path to excel metadata file',
        required=True)
    args = argp.parse_args()

    run_discrimination(args.meta, args.fasta)

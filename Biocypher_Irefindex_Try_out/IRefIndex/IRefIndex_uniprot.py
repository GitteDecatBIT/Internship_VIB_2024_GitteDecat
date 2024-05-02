import pypath.resources.urls as urls
import pypath.share.curl as curl
import re
from typing import Iterable

# Specify organism 
url = 'https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/7227.mitab.08-28-2023.txt.zip'
c = curl.Curl(url, silent = False, large = True, slow = True)
f = next(iter(c.result.values()))
nul = f.readline()
for l in f:
    l = l.split('\t')
    # ORGANISM
    input_organism= l[10]
    pattern_organism= r'taxid:(\d+)'
    match_organism= re.search(pattern_organism, input_organism)
    if match_organism:
        organism = match_organism.group(1) # Extracting the number
    else:
        organism = ""


class UniprotQuery:

    _PROCESS = {
        'dict': '_process_dict',
        'list': '_process_list',
    }
    _OP = ('_AND', '_NOT', '_OR')
    _OPSTART = re.compile(r'^(OR|AND)')
    _OPEND = re.compile(r'(OR|AND)$')
    _FIELDSEP = re.compile(r'[\s;]')
    _FIELDEND = re.compile(r'$;')
    _SYNONYMS = {
        'organism': 'organism_id',
        'ncbi_tax_id': 'organism_id',
    }
    _FIELD_SYNONYMS = {
        'function': 'cc_function',
        'activity_regulation': 'cc_activity_regulation',
        'tissue_specificity': 'cc_tissue_specificity',
        'developmental_stage': 'cc_developmental_stage',
        'induction': 'cc_induction',
        'intramembrane': 'ft_intramem',
        'signal_peptide': 'ft_signal',
        'subcellular_location': 'cc_subcellular_location',
        'transmembrane': 'ft_transmem',
        'comment': 'cc_miscellaneous',
        'topological_domain': 'ft_topo_dom',
        'family': 'protein_families',
        'interactor': 'cc_interaction',
        'keywords': 'keyword',
    }


def _swissprot_param(swissprot):

    return (
        'true'
            if swissprot in {'true', 'True', 'yes', 'YES', True} else
        'false'
            if swissprot in {'false', 'False', 'no', 'NO', False} else
        None
    )


def _all_uniprots(organism = organism, swissprot = None):

    swissprot = _swissprot_param(swissprot)
    rev = '' if swissprot is None else ' AND reviewed: %s' % swissprot
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism_id:%s%s' % (str(organism), rev),
        'format': 'tsv',
        'fields': 'accession',
    }

    if organism == '*':
        get['query'] = rev.strip(' AND ')

    c = curl.Curl(url, get = get, silent = False, slow = True)
    data = c.result

    return {
        l.strip() for l in data.split('\n')[1:] if l.strip()
    }

def uniprot_data(
        *query,
        fields: str | Iterable[str] | None = None,
        organism: str | int | None = 9606,
        reviewed: bool | None = True,
        **kwargs
    ) -> dict[str, str] | dict[str, dict[str, str]]:
    """
    Basic client for the UniProt REST API.

    Retrieves one or more fields from UniProt, by default for all reviewed
    (SwissProt) proteins of one organism

    Args:
        query:
            Query elements: can be a ready query or its components, bypassing
            the processing in this function or performing only simple
            concatenation. Alternatively, it can be a nested structure of lists
            and dicts describing more complex queries. See the examples below.
        fields:
            One or more UniProt field name. See details.
        organism:
            Organism name or identifier, e.g. "human", or "Homo sapiens",
            or 9606.
        reviewed:
            Restrict the query to SwissProt (True), to TrEMBL (False), or
            cover both (None).
        kwargs:
            Same as passing a dict to ``query``.

    Details:
        The query can be built in several ways:
        - Simple string or concatenation of strings:
          query_builder('kinase AND organism_id:9606')
          query_builder('kinase', 'organism_id:9606')
          query_builder('kinase', organism_id = 9606)
          The above 3 examples all return the same query:
          `kinase AND organism_id:9606`
        - The default operator within lists is `OR` and within dicts is `AND`:
          query_builder(organism = [9606, 10090, 10116])
          `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
          query_builder({'organism_id': 9606, 'reviewed': True})
          `organism_id:9606 AND reviewed:true`
        - These default operators can be changed by including the `op` key in
          dicts or including the operator with an underscore in lists:
          query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
          `length:[500 TO *] OR mass:[50000 TO *]`
          query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
          `lit_author:Huang AND lit_author:Kovac`
        - The nested structures translate into nested parentheses in the query:
          query_builder({'organism_id': [9606, 10090], 'reviewed': True})
          `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
        - Values are converted to strings, intervals can be provided as tuples:
          query_builder({'length': (100, None), 'organism_id': 9606})
          `length:[100 TO *] AND organism_id:9606`

        For a complete reference of the available parameters, see
        https://www.uniprot.org/help/query-fields and
        https://www.uniprot.org/help/text-search for additional syntax
        elements.

        For the available fields refer to the ``_FIELD_SYNONYMS`` attribute of
        the UniprotQuery class or the UniProt website:
        https://www.uniprot.org/help/return_fields

    Returns:
        - A list of UniProt IDs if no fields were provided.
        - A dict of UniProt IDs and corresponding field values if
          exactly one field was provided.
        - A dict with field names as top level keys and dicts of the
          kind described in the previous point as values.
    """

    for arg in ('organism', 'reviewed'):

        if locals()[arg] is not None:

            kwargs[arg] = locals()[arg]

    return uniprot_query(*query, fields = fields, **kwargs)


def uniprot_query(
        *query,
        fields: str | Iterable[str] | None = None,
        **kwargs
    ) -> dict[str, str] | dict[str, dict[str, str]]:
    """
    Basic client for the UniProt REST API.

    Args:
        query:
            Query elements: can be a ready query or its components, bypassing
            the processing in this function or performing only simple
            concatenation. Alternatively, it can be a nested structure of lists
            and dicts describing more complex queries. See the examples below.
        fields:
            One or more UniProt field name. See details.
        kwargs:
            Same as passing a dict to ``query``.

    Details:
        The query can be built in several ways:
        - Simple string or concatenation of strings:
          query_builder('kinase AND organism_id:9606')
          query_builder('kinase', 'organism_id:9606')
          query_builder('kinase', organism_id = 9606)
          The above 3 examples all return the same query:
          `kinase AND organism_id:9606`
        - The default operator within lists is `OR` and within dicts is `AND`:
          query_builder(organism = [9606, 10090, 10116])
          `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
          query_builder({'organism_id': 9606, 'reviewed': True})
          `organism_id:9606 AND reviewed:true`
        - These default operators can be changed by including the `op` key in
          dicts or including the operator with an underscore in lists:
          query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
          `length:[500 TO *] OR mass:[50000 TO *]`
          query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
          `lit_author:Huang AND lit_author:Kovac`
        - The nested structures translate into nested parentheses in the query:
          query_builder({'organism_id': [9606, 10090], 'reviewed': True})
          `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
        - Values are converted to strings, intervals can be provided as tuples:
          query_builder({'length': (100, None), 'organism_id': 9606})
          `length:[100 TO *] AND organism_id:9606`

        For a complete reference of the available parameters, see
        https://www.uniprot.org/help/query-fields and
        https://www.uniprot.org/help/text-search for additional syntax
        elements.

        For the available fields refer to the ``_FIELD_SYNONYMS`` attribute of
        the UniprotQuery class or the UniProt website:
        https://www.uniprot.org/help/return_fields

    Returns:
        - A list of UniProt IDs if no fields were provided.
        - A dict of UniProt IDs and corresponding field values if
          exactly one field was provided.
        - A dict with field names as top level keys and dicts of the
          kind described in the previous point as values.
    """

    return UniprotQuery(*query, fields = fields, **kwargs).perform()


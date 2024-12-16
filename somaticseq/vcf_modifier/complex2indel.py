#!/usr/bin/env python3


def translate(refbase, altbase):
    offset = 0

    if len(refbase) == len(altbase):
        return False

    elif len(refbase) == 1 or len(altbase) == 1:
        return ((refbase, altbase), offset)

    else:
        for base_i, base_j in zip(refbase[::-1], altbase[::-1]):
            if base_i == base_j and (len(refbase) >= 2 and len(altbase) >= 2):
                refbase = refbase[:-1]
                altbase = altbase[:-1]
            else:
                break

        for base_i, base_j in zip(refbase, altbase):
            if base_i == base_j and (len(refbase) >= 2 and len(altbase) >= 2):
                refbase = refbase[1:]
                altbase = altbase[1:]
                offset += 1
            else:
                break

        return ((refbase, altbase), offset)


def resolve_complex_variants_into_snvs_and_indels(
    refbases: str, altbases: str
) -> list[dict] | None:
    """
    Split complex variants into combination of snvs and indels.
    """

    if len(refbases) == 1 or len(altbases) == 1:  # snv / indel
        return None

    # Initialize a list to hold the new records
    list_of_variants: list[dict] = []

    # "Left-align" the REF and ALT to assign snvs until one has to consider
    # deletion or insertion
    for i, (refbase, altbase) in enumerate(zip(refbases, altbases)):
        if refbase != altbase:
            list_of_variants.append(
                {"OFFSET": i, "REF": refbase, "ALT": altbase},
            )
    # Handle deletion
    if len(refbases) > len(altbases):
        list_of_variants.append(
            {"OFFSET": i, "REF": refbases[i:], "ALT": altbases[i]},
        )
    # Handle insertion
    elif len(altbases) > len(refbases):
        list_of_variants.append(
            {"OFFSET": i, "REF": refbases[i], "ALT": altbases[i:]},
        )
    return list_of_variants

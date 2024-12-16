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

def resolve_complex_variants_into_snvs_and_indels(
    refbases: str, altbases: str
) -> list[dict]:
    """
    Split complex variants into combination of snvs and indels.
    """
    snv_or_indel = [{"OFFSET": 0, "REF": refbases, "ALT": altbases}]

    if len(refbases) == 1 and len(altbases) == 1:  # snv
        return snv_or_indel

    if (len(refbases) == 1 or len(altbases) == 1) and (
        refbases[0] == altbases[0]
    ):  # indel
        return snv_or_indel

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
            {"OFFSET": i, "REF": refbases[i:], "ALT": refbases[i]},
        )
    # Handle insertion
    elif len(altbases) > len(refbases):
        list_of_variants.append(
            {"OFFSET": i, "REF": refbases[i], "ALT": refbases[i] + altbases[i + 1 :]},
        )
    return list_of_variants

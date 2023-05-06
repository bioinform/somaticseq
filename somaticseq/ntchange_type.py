def ntchange(variant_frame):
    GC2CG = []
    GC2TA = []
    GC2AT = []
    TA2AT = []
    TA2GC = []
    TA2CG = []

    for ref, alt in zip(variant_frame["REF"], variant_frame["ALT"]):
        ref = ref.upper()
        alt = alt.upper()

        if (ref == "G" and alt == "C") or (ref == "C" and alt == "G"):
            GC2CG.append(1)
            GC2TA.append(0)
            GC2AT.append(0)
            TA2AT.append(0)
            TA2GC.append(0)
            TA2CG.append(0)

        elif (ref == "G" and alt == "T") or (ref == "C" and alt == "A"):
            GC2CG.append(0)
            GC2TA.append(1)
            GC2AT.append(0)
            TA2AT.append(0)
            TA2GC.append(0)
            TA2CG.append(0)

        elif (ref == "G" and alt == "A") or (ref == "C" and alt == "T"):
            GC2CG.append(0)
            GC2TA.append(0)
            GC2AT.append(1)
            TA2AT.append(0)
            TA2GC.append(0)
            TA2CG.append(0)

        elif (ref == "T" and alt == "A") or (ref == "A" and alt == "T"):
            GC2CG.append(0)
            GC2TA.append(0)
            GC2AT.append(0)
            TA2AT.append(1)
            TA2GC.append(0)
            TA2CG.append(0)

        elif (ref == "T" and alt == "G") or (ref == "A" and alt == "C"):
            GC2CG.append(0)
            GC2TA.append(0)
            GC2AT.append(0)
            TA2AT.append(0)
            TA2GC.append(1)
            TA2CG.append(0)

        elif (ref == "T" and alt == "C") or (ref == "A" and alt == "G"):
            GC2CG.append(0)
            GC2TA.append(0)
            GC2AT.append(0)
            TA2AT.append(0)
            TA2GC.append(0)
            TA2CG.append(1)

        else:
            GC2CG.append(0)
            GC2TA.append(0)
            GC2AT.append(0)
            TA2AT.append(0)
            TA2GC.append(0)
            TA2CG.append(0)

    new_data = variant_frame.assign(
        GC2CG=GC2CG, GC2TA=GC2CG, GC2AT=GC2CG, TA2AT=GC2CG, TA2GC=GC2CG, TA2CG=GC2CG
    )
    return new_data

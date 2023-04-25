#!/usr/bin/env python3


class BedFile:
    def __init__(self, BedFile):

        """Argument is a line in pileup file."""
        self.BedFile = BedFile

        bedRegions = {}

        with open(self.BedFile) as bed:
            line_i = bed.readline().rstrip()

            while line_i:
                item = line_i.split("\t")

                contig = item[0]
                region = (int(item[1]), int(item[2]))

                if contig not in bedRegions:
                    bedRegions[contig] = []

                bedRegions[contig].append(region)

                line_i = bed.readline().rstrip()

        self.bedRegions = bedRegions

    def inRegion(self, contig_i, position_i, ordered=True):

        intersected = False

        if contig_i in self.bedRegions:

            # If the BED file is ordered, it can break out once it goes beyond the position_i
            if ordered:
                for region_i in self.bedRegions[contig_i]:

                    if region_i[0] < position_i <= region_i[1]:
                        intersected = True
                        break

                    elif region_i[0] > position_i:
                        break

            # If the BED file is not ordered, then it needs to go all the way to the end every time
            else:
                for region_i in self.bedRegions[contig_i]:

                    if region_i[0] < position_i <= region_i[1]:
                        intersected = True
                        break

        return intersected

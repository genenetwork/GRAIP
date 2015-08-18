"""
This application generates permuted genomes for QTL analysis of an AIL.

Instead of permuting phenotypes, this approach to the problem of permuting
an AIL permutes the identity of the next to last generation (the parents
of the final generation) and regenerates the genotype file for a final
generation population based on these permuted identities.

File formats:
Map File
One line for each chromosome, with marker positions in cM

Haplotypes File
The following slightly unusual format:
Each parent is a column.  The first row is name.  Second row is sex
Third and subsequent rows are genotypes, with the entire set of genotypes
for strand one first, then all genotypes for strand 2.
Assumption is that X is the last chromosome

Family
Each "family" (sibship) in the family file is a separate line with
momID, dadID, number_of_offspring separated by white space
followed by the sexes of the offspring separated by white space.
The easiest way to specify offspring is 1 per line, but many per line
is acceptable.
"""


from random import shuffle, sample, uniform
"""
Shuffle randomly re-orders a list.
Sample samples without replacement from a list.
Uniform samples from a uniform distribution of range specified.
"""

from RandomArray import poisson
"""poisson(mean) draws items from a poisson distribution with an average equal
to mean.  This function is part of Numeric23.8"""

class AILPerm(object):
    # Codes for various genotype combinations in the output file.
    _genoComboLookup = {("B","B"):"0", ("D","D"):"2", ("B","D"):"1", \
    ("D","B"):"1", ("B","U"):"9",("U","B"):"9", ("D","U"):"9", ("U","D"):"9", \
    ("U","U"):"9"}
    
    def loadMap(self, filename):
        """Uses map file to populate self.map and self.chrLength

        filename is a string that is the name of the file holding the map data
        in the format described above."""
        self.map = []
        self.chrLength = []

        """gets the lines in the file and changes their elements to numbers."""
        lines = splitFileIntoLists(filename)
        for line in lines:
            self.map.append(map(float, line))

        """ here we also find lengths of each chromosome (distance to last
        marker). If there are chiasmata before the first marker that will only
        flip the "baseStrand" which doesn't matter b/c it is arbitrary to
        begin with."""
        for chr in self.map:
            self.chrLength.append(max(chr))


    def loadParents(self, filename):
        """loads in haplotypes from haplotypes file described above"""
        self.dads = [] #initializes the list we will use for the dads
        self.moms = [] #initializes the list we will use for the moms
        markersByChr = [0]
        total_markers = 0
        
        lines = splitFileIntoLists(filename)    #reads file.
        parents = map(list, zip(*lines))
        """each parent is now  a list each element of which has
        parent_id, sex, strand elements, other strand elements"""
        strand_length = (len(parents[0])- 2)/2

        """this loop makes a list of 1 number per chr which are the cumulative
        number of markers for each chr, used later"""         
        for chr in self.map:
            total_markers += len(chr) 
            markersByChr.append(total_markers)

        """This loop runs through each parent list (from above) and separates
        into strands, then adds that to either the moms or dads file"""        
        for parent in parents:
            parentGenos = parent[2:] #snips off parentID and sex
            strandOne = []
            strandTwo = []
            expStrandLen = 2*strand_length + 2

            if len(parent)<> expStrandLen: #check to see if strand length is OK
                raise ValueError, "Error in length of genotype set.  \
                Everything else is innacurate"

            """For the number of markers in each chromosome append a list for
            each chromosome to each strand.  At this point strandOne and
            strandTwo are lists of lists of genotypes by chromosome."""
            for markerNum in markersByChr:
                if markerNum > 0: #the first time through the loop, just update last_markerNum.  The next times do the loop
                    strandOne.append(parentGenos[last_markerNum:markerNum]) #cuts up parent list for strand 1 into chromosomes
                    strandTwo.append(parentGenos[last_markerNum+strand_length:markerNum+strand_length])
                last_markerNum = markerNum

            """add to dads or moms file as appropriate"""
            if parent[1] == "m" or parent[1]=='M':  #makes file of dads.  each dad is a list of the form [dad_number, [strand 1], [strand 2]]
                self.dads.append([parent[0],[strandOne, strandTwo]])
            if parent[1] == "f" or parent[1]=='F':  #makes file of dads as per above
                self.moms.append([parent[0],[strandOne, strandTwo]])


    def loadFamilies(self, filename):
        """Loads each family in the proper format: momID, dadID, n_offspr, list of sexes of offspring.
        One nice hint is that if you list each offspring on a separate line (so n_offspr=1) this
        becomes a nice rectangular table."""
        
        families = []
        familyFile = splitFileIntoLists(filename)

        """Formats families into a list of lists, with each family being
        in the format [momID, dadID, number_of_offspring, sexes of offspring]"""
        for family in familyFile:
            formattedFamily = [family[0], family[1], family [2], family[3:]]
            families.append(formattedFamily)
        return families 


    def permuteParents(self):
        """this function permutes the parent IDs and creates a
        dictionary entry for each."""
        
        self.lookupMom = {}
        self.lookupDad = {}

        """zips, shuffles, and returns in the original format the list of dads,
        permuting the parental ID labels"""
        dads = zipShuffleReturn(self.dads)
        for dad in dads:
            self.lookupDad[dad[0]] = dad[1]

        """zips, shuffles, and returns in the original format the list of dads,
        permuting the parental ID labels"""    
        moms = zipShuffleReturn(self.moms)
        for mom in moms:
            self.lookupMom[mom[0]] = mom[1]


    def recombine(self, parentID, sex):
        """Given a parentID and the sex of the parent, this function looks
        up the haplotype information for that parent and generates from
        that parent a recombined haplotype that is inherited by a zygote."""

        """this is where the alleles at each marker will be placed, as a
        simple list of alleles."""
        haplotype = []
        
        if sex == 'm':
            parent = self.lookupDad[parentID]
        if sex == 'f':
            parent = self.lookupMom[parentID]

        """makes a set of chiasmata that are overlaid on the haplotype
        causing switching of strands at each breakpoint"""
        chiasmata = self.makeChiasmata()
        
        numChr = len(parent[0]) # numChr is the number of chromosomes

        """for each chromosome 0 to chrNum this loop returns a haplotype with
        overlaid chiasmata. Strands are in parent[0] and parent[1]."""
        for chrNum in xrange(numChr):
            baseStrandID = sample((0,1),1)[0] #randomly choose a base strand
            alleleNum = 0
            if baseStrandID == 0:
                baseStrand = parent[:] #strand order is the default
            if baseStrandID == 1:
                baseStrand = [parent[1],parent[0]] #strand order is reversed
                
            """if no chiasmata, reads the contents of the chromosome into the haplotype."""
            if chiasmata[chrNum] == []:
                for allele in baseStrand[0][chrNum]:
                    haplotype.append(allele)

            """if there are chiasmata, flip the strand you're reading from at
            each one"""
            else:
                for marker in self.map[chrNum]:
                    
                    """counter for number of chiasmata up to that point"""
                    previousChiasmata = 0

                    """this loop counts the number of breakpoints previous to
                    the marker"""
                    for breakpoint in chiasmata[n]: 
                        if breakpoint<marker:
                            previousChiasmata +=1

                    """We use the modulus 2 value -- 0 or 1 for
                    even/odd numbers of chiasmata to chose the strain based on
                    how many chiasmata have passed."""
                    strandNum = previousChiasmata%2  #  %2 means mod 2

                    """appends the correct haplotype strand to haplotype"""
                    haplotype.append(baseStrand[strandNum][chrNum][alleleNum])

                    """alleleNum is iterated just before marker to update the
                    position of the allele that is inserted"""
                    alleleNum += 1
        return haplotype

                        
    def makeChiasmata(self):
        """makes a list of chiasmata for each chromosome using the simple model outlined below"""
        
        interference = 5 # absolute interference distance in cM units

        """number of times to check if chiasmata are too close to each other.
        Since there are only #a few chiasmata per chromosome this doesn't need
        to be done all that many times."""
        draws = 100 
                    
        genomeChiasmata = []
        """this is a very simple model.  it uses a poisson model to determine
        how many recombinations are seen on a chromosome a with a total length
        of N cM If there are any, it draws that number of positions from a
        uniform distribution. If there are more than one, it sorts them by
        position, and checks to see whether they are close to each other by
        choosing one randomly. If they are within the (absolute) interference
        distance one of the two (randomly chosen) is discarded.
        
        This process is repeated enough times to ensure that all chiasma are
        at least that radius from each other."""
        
        for chr in self.chrLength: # chr starts in cM units
            chromosomeChiasmata = []
            
            #devines how many chiasmata there will be given the chr length.
            poissonEvents = poisson(chr/100)

            if poissonEvents == 0:
                genomeChiasmata.append([])

            else:
                """make a list of uniform draws of poisson events"""
                for event in xrange(poissonEvents):  
                    chromosomeChiasmata.append(uniform(1,chr))
                chromosomeChiasmata.sort() #sort the list

                for n in xrange(draws): # make pairs and compare them.
                    """if there is only one, break out of the loop."""
                    if len(chromosomeChiasmata)==1:
                        break

                    """for more than one snip out one at random if it closer to
                    the other than is allowed by the interference variable"""
                    index = sample(range(len(chromosomeChiasmata)),2)
                    #if closer than interference param
                    if abs(chromosomeChiasmata[index[0]]-chromosomeChiasmata[index[1]])<interference:
                        #snip out one
                        chromosomeChiasmata = chromosomeChiasmata[:index[0]]+chromosomeChiasmata[index[0]+1:]                 genomeChiasmata.append(chromosomeChiasmata)
        return genomeChiasmata


    def makeZygote(self, egg, sperm, sex):
        """takes the recombined (or not) haplotypes and puts them together as a
        zygote using the genotype coding from genoComboLookup"""

        """num markers on X, assuming X is the last chromosome"""
        numXMarkers = len(self.map[len(self.map)-1]) #length of final chr.
        
        lastAutosomalMarker = len(egg) - 1 - numXMarkers

        doubleStrand = []
        zygote = []

        """Handles appending egg and sperm genotypes.  If the animal is a
        male and the chromosome is X, appends two copies of the egg X
        genotype."""
        for n in xrange(len(egg)):
            if sex == 'M' and n > lastAutosomalMarker: #if it's the X
                doubleStrand.append((egg[n],egg[n]))
            else:
                doubleStrand.append((egg[n],sperm[n])) #if it's not the X

        """looks up each pair of egg, sperm genotypes and returns zygote geno"""        
        for pair in doubleStrand:
            zygote.append(AILPerm._genoComboLookup[pair])
            
        return zygote


    def makePop(self, mapFile, genoFile, familyFile):
        """control function.  Loads the files, runs the permutation, then for each permuted family
        it makes n offspring."""

        """loads files"""
        self.loadMap(mapFile)
        self.loadParents(genoFile)
        families = self.loadFamilies(familyFile)

        """permutes parents"""
        self.permuteParents()
        
        population = []

        """For each family, this gets the mom and dad genotypes, directs
        generation of sperm and egg haplotypes, and zygote genotype
        and appends each zygote to the population"""
        for family in families:
            momID = family[0]
            dadID = family[1]
            n_offspring = int(family[2]) #number of offspring
            sexList = family[3] # list of sexes of offspring, often a single M/F

            """for each offspring, make an egg and sperm, then a zygote, and
            add to the population"""
            for i in xrange(n_offspring):
                egg = self.recombine(momID,'f')
                sperm = self.recombine(dadID,'m')
                zygote = self.makeZygote(egg, sperm, sexList[i])
                population.append(zygote)
                
        return population


def makeOutput(input, fileName):
    """turns list of lists into a \n delim string"""
    listForJoining = []
    for line in input:
        listForJoining.append(" ".join(line))
    output = "\n".join(listForJoining) + "\n"
    file = open(fileName,mode='w')
    file.write(output)
    file.close()


def splitFileIntoLists(filename):
    """reads in a file, splits each line on white space"""
    rows = []
    a_file = open(filename)
    for line in a_file:           #each line is an entry
        rows.append(line.split())
    return rows


def zipShuffleReturn(listname):
    """zips first elements of a list together, shuffles, rezips, and returns"""
    listToReturn = map(list, zip(*listname))
    shuffle(listToReturn[0])
    listToReturn = map(list, zip(*listToReturn))
    return listToReturn


if __name__ == "__main__":
    try:
        myAILPerm = AILPerm()
        start = 1                 #the program will generate files *including* the start and endpoint.  start=1, stop=2
        stop = 2                  #for instance will generate geno1.txt, geno2.
        
        """These are the data files described at the top"""
        mapData = 'ailmap.txt'
        haploData = 'ailhaplo.txt'
        famData = 'ailfam.txt'

        """Makes one population for each number, inclusive, between start and stop"""
        for n in xrange(start,stop+1):
            population = myAILPerm.makePop(mapData,haploData,famData)
            makeOutput(population, "geno"+str(n)+".txt")

    except Exception, e:
        print 'Program encountered an error', e







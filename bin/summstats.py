import sys, os, egglib, random, math

class Application(object):
    """
    Class representing the program. The class is used in a very simple
    manner::

        >>> app = Application()
        >>> app.run(args)

    Where ``args`` is the list of arguments (usually: ``sys.args[1:]``).
    """

    manual = """summstats

1. Command line arguments

The program is a single Python script. It is invoked as this:
    $ python summstats.py <OPTIONS>

Options are all in the form <KEY>=<VALUE> with keys as listed below:
    - input-file        name of the input file
    - output-file       name of the output file
    - LSS               comma-separated list of summary statistics
                        computed on selected loci, taken from the list:
                            * He          heterozygosity
                            * Dj          Jost's D
                            * WCst        Weir and Cockerham's theta
    - WSS               list of statistics computed on a window around
                        selected loci, taken LSS plus that list:
                            * S           number of polymorphic sites
                            * thetaW      Watterson's 4Nu estimator
                            * D           Tajima's D
                            * Da          net distance between populations
                            * ZnS         Kelly et al.'s ZnS
                            * ZZ          Rozas et al.'s ZZ
    - GSS               list of global summary statistics, taken from
                        the WSS list, plus SFS
    - GSSp              global summary statistics that should be computed
                        for each population (excluding SFS, Dj, and WCst)
    - wspan             span of windows on each side of the focal locus
                        (in bp, such that the window size is actually
                        2 * wspan + 1 bp)
    - SFS-bins          number of bins for SFS
    - select            method to select focal loci (to compute LSS and
                        WSS), from the list:
                            * all         use all loci
                            * rand        use each locus randomly
                            * freq        use each nth locus
                            * list        use the "selection" column of
                                          the input file
    - select-num        number of loci to draw randomly (if select=rand)
    - select-freq       period of selected loci (if select=freq): if
                        freq=1, all loci are selected; if freq=2, each
                        second locus is selected; and so on
    - debug             activate debug mode:
                            * 0   default mode
                            * 1   debug mode: exception trace displayed
                        when setting this option to 1, it is recommended
                        to do it as first argument to cope with an error
                        occurring during argument parsing 

2. Input file

The input file is a text file with space/tab separated values. There is
a header line giving column names and, then, one line per site (a site
may be a SNP, and indel or whatever). No term may contain a space. The
order of columns is unimportant, except that the individual columns
must be after the columns with fixed names.

List of columns:
    "chrom"       a string identifying a chromosome, contig or whatever.
    "position"    site position, as an integer.
    "status"      site category, chosen from "IG" (intergenic), "S"
                  (synonymous), and "NS" (non-synonymous).
    "alleles"     a comma-separated list of allelic values. Each allele
                  is represented by a string or undefined length. There
                  must be at least two alleles.
    "selection"   a column of Y or N to indicate whether LSS statistics
                  should be reported for each site [OPTIONAL].
    "name@label"  one additional column per individual, where "name" is
                  the sample name and "label" is a population label
                  (both are strings).
    ...           and so on for all individuals.

The data (values for individual columns) are two-digit strings with
integer values (diploid data only). 
    - 0 for missing data.
    - 1 for the first allele.
    - 2 for the second allele, and so on.

3. Output file

The output file is a tab-separated table, with a header line and, then,
one line per locus (matching the loci found in the input file.

The first column is titled "ID" and gives the identifier of a loci, in
the form "chrom:position", as taken from the input file.

All subsequent columns provide the values for a statistic. The column
name is the statistic code (with a "LSS", "WSS", or "GSS" prefix). For
GSS, the value of the statistics is identical for all loci. Missing data
are denoted by "NA".
"""

    valid_LSS = ['He', 'Dj', 'WCst']
    valid_LSSp = ['He', 'Dj']
    valid_WSS = valid_LSS + ['S', 'thetaW', 'Pi', 'D', 'Da', 'ZZ', 'ZnS']
    valid_WSSp = valid_LSSp + ['S', 'thetaW', 'Pi', 'D', 'ZZ', 'ZnS']
    valid_GSS = valid_WSS[:]
    valid_GSSp = ['He', 'S', 'thetaW', 'Pi', 'D', 'Da', 'ZZ', 'ZnS']
    valid_select = ['all', 'rand', 'freq', 'list']
    alph = egglib.alphabets.Alphabet('int', [1,2,3], [-1])

    def run(self, args):
        """
        Run the program, using the provided list of arguments.
        """

        if len(args) == 0:
            print self.manual
        else:
            self.get_args(args)
            self.get_input_file()
            self.compute_SFS()
            self.compute_GSS()
            self.compute_GSSp()
            self.compute_LSS()
            self.compute_LSSp()
            self.compute_WSS()
            self.compute_WSSp()
            self.export()

    def get_args(self, args):
        """
        This method interprets the list of arguments, check their
        consistency, and set the ``_params`` dictionary. In particular,
        lists of statistics are in ``_params['GSS']``,
        ``_params['LSS']``, and ``_params['WSS']`` (as well as their
        per-population counterparts).
        """

        # initialize parameter dict
        self._params = {'LSS': set(), 'LSSp': set(), 'WSS': set(),
                'WSSp': set(), 'GSS': set(), 'GSSp': set(),
                'SFS': False, 'debug': 0}

        # get arguments
        for arg in args:
            if arg.count('=') != 1:
                raise ValueError, 'invalid command-line argument: \'{0}\''.format(arg)
            key, value = arg.split('=')
            if key == 'input-file':
                self._params['input-file'] = value
            elif key == 'output-file':
                self._params['output-file'] = value
            elif key == 'LSS':
                stats = value.split(',')
                for stat in stats:
                    if stat not in self.valid_LSS: raise ValueError, 'invalid LSS: {0}'.format(stat)
                self._params['LSS'].update(stats)
            elif key == 'LSSp':
                stats = value.split(',')
                for stat in stats:
                    if stat not in self.valid_LSSp: raise ValueError, 'invalid LSSp: {0}'.format(stat)
                self._params['LSSp'].update(stats)
            elif key == 'GSS':
                stats = value.split(',')
                for stat in stats:
                    if stat == 'SFS':
                        self._params['SFS'] = True
                    elif stat not in self.valid_GSS: raise ValueError, 'invalid GSS: {0}'.format(stat)
                self._params['GSS'].update(stats) # SFS is added but will be filtered out during sorting
            elif key == 'GSSp':
                stats = value.split(',')
                for stat in stats:
                    if stat not in self.valid_GSSp: raise ValueError, 'invalid GSSp: {0}'.format(stat)
                self._params['GSSp'].update(stats)
            elif key == 'WSS':
                stats = value.split(',')
                for stat in stats:
                    if stat not in self.valid_WSS: raise ValueError, 'invalid WSS: {0}'.format(stat)
                self._params['WSS'].update(stats)
            elif key == 'WSSp':
                stats = value.split(',')
                for stat in stats:
                    if stat not in self.valid_WSSp: raise ValueError, 'invalid WSSp: {0}'.format(stat)
                self._params['WSSp'].update(stats)
            elif key == 'wspan':
                try: value = int(value)
                except ValueError: raise ValueError, 'invalid value for \'window-size\': {0}'.format(value)
                if value < 0: raise ValueError, '\'window-size\' must be >= 0'
                self._params['wspan'] = value
            elif key == 'select':
                if value not in self.valid_select: raise ValueError, 'invalud value for \'select\': {0}'.format(value)
                self._params['select'] = value
            elif key == 'select-num':
                try: value = int(value)
                except ValueError: raise ValueError, 'invalid value for \'select-num\': {0}'.format(value)
                if value < 0: raise ValueError, '\'select-num\' must be >=0'
                self._params['select-num'] = value
            elif key == 'select-freq':
                try: value = int(value)
                except ValueError: raise ValueError, 'invalid value for \'select-freq\': {0}'.format(value)
                if value < 1: raise ValueError, '\'select-freq\' must be >0'
                self._params['select-freq'] = value
            elif key == 'SFS-bins':
                try: value = int(value)
                except ValueError: raise ValueError, 'invalid value for \'SFS-bins\': {0}'.format(value)
                if value < 1: raise ValueError, '\'SFS-bins\' must be >0'
                self._params['SFS-bins'] = value
            elif key == 'debug':
                if value == '0': self._params['debug'] = 0
                elif value == '1': self._params['debug'] = 1
                else: raise ValueError, 'invalud value for \'debug\': {0}'.format(value)
            else:
                raise ValueError, 'invalid option: ' + key

        # check required parameters & consistency
        if 'input-file' not in self._params: raise ValueError, '\'input-file\' argument is required'
        if 'output-file' not in self._params: raise ValueError, '\'output-file\' argument is required'
        if 'select' not in self._params: raise ValueError, '\'select\' argument is required'
        if len(self._params['WSS']) > 0 and 'wspan' not in self._params: raise ValueError, '\'wspan\' argument is required because of required WSS'
        if len(self._params['WSS']) == 0 and 'wspan' in self._params: sys.stderr.write('warning: \'wspan\' is specified but won\'t be used\n')
        if self._params['select'] == 'rand' and 'select-num' not in self._params: raise ValueError, '\'select-num\' argument is required because \'select\' is \'rand\''
        if self._params['select'] != 'rand' and 'select-num' in self._params: sys.stderr.write('warning: \'select-num\' is specified but won\'t be used\n')
        if self._params['select'] == 'freq' and 'select-freq' not in self._params: raise ValueError, '\'select-freq\' argument is required because \'select\' is \'freq\''
        if self._params['select'] != 'freq' and 'select-freq' in self._params: sys.stderr.write('warning: \'select-freq\' is specified but won\'t be used\n')
        if self._params['SFS'] == True and 'SFS-bins' not in self._params: raise ValueError, '\'SFS-bins\' argument is required because \'SFS\' is requested'
        if self._params['SFS'] == False and 'SFS-bins' in self._params: sys.stderr.write('warning: \'SFS-bins\' is specified but won\'t be used\n')

        # order statistics
        self._params['GSS'] = [i for i in self.valid_GSS if i in self._params['GSS']]
        self._params['LSS'] = [i for i in self.valid_LSS if i in self._params['LSS']]
        self._params['WSS'] = [i for i in self.valid_WSS if i in self._params['WSS']]
        self._params['GSSp'] = [i for i in self.valid_GSSp if i in self._params['GSSp']]
        self._params['LSSp'] = [i for i in self.valid_LSSp if i in self._params['LSSp']]
        self._params['WSSp'] = [i for i in self.valid_WSSp if i in self._params['WSSp']]

    def get_input_file(self):
        """
        This methods reads the input file, checks that column
        `selection` is present if needed (and warns if it is present if
        not needed) and set the _sites, _positions, and _selection
        members.

        Store data in the following members: ``_id`` (list of site ID),
        ``_sites`` (list of egglib.stats.Site instances), ``_selection``
        (list of indexes of loci selected for LSS and WSS), ``_status``
        (list of site statuses), ``_pos`` (list of site positions),
        ``_chrom`` (list of chromosome for each site) and ``_struct``
        (an egglib.stats.Structure instance).
        """

        # open file and process header
        try: f = open(self._params['input-file'])
        except IOError: raise ValueError, 'cannot open this file: {0}'.format(self._params['input-file'])
        header = f.readline().split()
        if set(header[:4]) == set(['chrom', 'position', 'status', 'alleles']):
            selectionQ = False
            if self._params['select'] == 'list': raise ValueError, '\'selection\' column is required in input file because \'select\' is \'list\''
            offset = 4
        elif set(header[:5]) == set(['chrom', 'position', 'status', 'alleles', 'selection']):
            selectionQ = True
            if self._params['select'] != 'list': sys.stderr.write('warning: \'selection\' column is present in input file but won\'t be used\n')
            offset = 5
        else:
            raise IOError, 'invalid header line in file {0}'.format(self._params['input-file'])
        self._sample = {}
        for idx, idv in enumerate(header[offset:]):
            if idv.count('@') != 1: raise IOError, 'invalid sample specification: {0}'.format(idv)
            name, pop = idv.split('@')
            if pop not in self._sample: self._sample[pop] = []
            self._sample[pop].append(name)
        self._pops = sorted(self._sample)

        # report sample structure
        print '[summstats] sample structure'
        print 'number of populations: {0}'.format(len(self._sample))
        print 'number of samples: {0}'.format(sum(map(len, self._sample.itervalues())))
        for pop in sorted(self._sample):
            print '    {0}:\t{1}'.format(pop, len(self._sample[pop]))

        # generate structure objects
        acc = 0
        ingr = {}
        for i, p in enumerate(self._pops):
            ingr[i] = {}
            for j in xrange(len(self._sample[p])):
                ingr[i][2*(j+acc)] =(2*(j+acc),)
                ingr[i][2*(j+acc)+1] = (2*(j+acc)+1,)
            acc += len(self._sample[p])
        self._struct = egglib.stats.make_structure({0: ingr}, {})
        self._struct_per_pop = {}
        for k in self._sample:
            idx = self._pops.index(k)
            self._struct_per_pop[k] = egglib.stats.make_structure({0: {idx: ingr[idx]}}, {})

        # initialize tables
        if selectionQ: self._selection = []
        self._id = []
        self._chrom = []
        self._pos = []
        self._sites = []
        self._status = []
        cur_chr = None
        cur_pos = None

        # process all sites
        for idx, line in enumerate(f):
            bits = line.split()
            if len(bits) != len(header): raise IOError, 'invalid number of items for site number {0}'.format(idx+1)
            line = dict(zip(header, bits))

            # get position
            try: pos = int(line['position']) - 1
            except ValueError: raise IOError, 'invalid position for site number {0}'.format(idx)
            if pos < 0: raise IOError, 'invalid position for site number {0}'.format(idx)
            if line['chrom'] != cur_chr:
                cur_chr = line['chrom']
                cur_pos = pos
            elif pos <= cur_pos:
                raise IOError, 'invalid ordering of sites detected at site number {0}'.format(idx)

            # set ID
            self._id.append('{0}:{1}'.format(cur_chr, pos+1))
            self._chrom.append(cur_chr)
            self._pos.append(pos)

            # add site to selection list
            if selectionQ:
                if line['selection'] not in set('YN'): raise IOError, 'invalid `selection` value at site number {0}'.format(idx)
                if line['selection'] == 'Y':
                    self._selection.append(len(self._id)-1)

            # get status
            if line['status'] not in set(['IG', 'NS', 'S']): raise IOError, 'invalid `status` value at site number {0}'.format(idx)
            self._status.append(line['status'])

            site = map(self.get_snp, [line['{0}@{1}'.format(sam, pop)] for pop in self._pops for sam in self._sample[pop]])
            site = [j for i in site for j in i]
            site = egglib.stats.site_from_list(site, [], self.alph)
            self._sites.append(site)

        # define list of sites
        if self._params['select'] == 'rand':
            self._selection = random.sample(xrange(len(self._sites)), self._params['select-num'])
            self._selection.sort()
        elif self._params['select'] == 'all':
            self._selection = xrange(len(self._sites))
        elif self._params['select'] == 'freq':
            self._selection = xrange(0, len(self._sites), self._params['select-freq'])

    def get_snp(self, x):
        """
        Convert genotype
        """
        if x == '0': return (-1, -1)
        if x == '00': return (-1, -1)
        if len(x) != 2: raise ValueError, 'invalid genotype: {0}'.format(x)
        try: a, b = map(int, x)
        except ValueError: raise ValueError, 'invalid genotype: {0}'.format(x)
        if a < 1 or b < 1: raise ValueError, 'invalid genotype: {0}'.format(x)
        return (a, b)

    def compute_SFS(self):
        """
        Compute the SFS. Store data in ``_SFS``.
        """
        if self._params['SFS']:
            self._SFS = [0] * self._params['SFS-bins']
            cs = egglib.stats.ComputeStats(multi_hits=True)
            cs.add_stats('MAF')
            for site in self._sites:
                MAF = cs.process_site(site)['MAF']
                if MAF is not None:
                    idx = int(math.ceil(MAF/0.5 * self._params['SFS-bins']))-1
                    self._SFS[idx] += 1

    def compute_GSS(self):
        """
        Compute global summary statistics. Store data in ``_GSS_stats``.
        """
        if len(self._params['GSS']) == 0:
            self._GSS_stats = []
            return
        cs = egglib.stats.ComputeStats(multi_hits=True, struct=self._struct)
        cs.add_stats(*self._params['GSS'])
        self._GSS_stats = cs.process_sites(self._sites)

    def compute_GSSp(self):
        """
        Compute global summary statistics per population. Store data in ``_GSSp_stats``.
        """
        self._GSSp_stats = {}
        if len(self._params['GSSp']) == 0: return
        cs = egglib.stats.ComputeStats()
        cs.add_stats(*self._params['GSSp'])
        for pop in self._pops:
            cs.configure(multi_hits=True, struct=self._struct_per_pop[pop])
            self._GSSp_stats[pop] = cs.process_sites(self._sites)

    def compute_LSS(self):
        """
        Compute locus summary statistics. Store data in ``_LSS_stats``.
        """
        self._LSS_stats = []
        if len(self._params['LSS']) == 0: return
        cs = egglib.stats.ComputeStats(multi_hits=True, struct=self._struct)
        cs.add_stats(*self._params['LSS'])
        for idx in self._selection:
            self._LSS_stats.append(cs.process_site(self._sites[idx]))

    def compute_LSSp(self):
        """
        Compute locus summary statistics. Store data in ``_LSSp_stats``.
        """
        self._LSSp_stats = {}
        if len(self._params['LSSp']) == 0: return
        cs = egglib.stats.ComputeStats()
        cs.add_stats(*self._params['LSSp'])
        for k in self._struct_per_pop: self._LSSp_stats[k] = []
        for idx in self._selection:
            for pop in self._pops:
                cs.configure(multi_hits=True, struct=self._struct_per_pop[pop])
                self._LSSp_stats[pop].append(cs.process_site(self._sites[idx]))

    def compute_WSS(self):
        """
        Compute window summary statistics. Store data in ``_WSS_stats``.
        """
        self._WSS_stats = []
        if len(self._params['WSS']) == 0: return
        cs = egglib.stats.ComputeStats(multi_hits=True, struct=self._struct)
        cs.add_stats(*self._params['WSS'])
        for idx in self._selection:
            first = idx
            while first > 0:
                if self._chrom[first-1] != self._chrom[idx] or self._pos[first-1] < self._pos[idx]-self._params['wspan']:
                    break
                first -= 1
            last = idx
            while last < len(self._sites) - 1:
                if self._chrom[last+1] != self._chrom[idx] or self._pos[last+1] > self._pos[idx]+self._params['wspan']:
                    break
                last += 1
            self._WSS_stats.append(cs.process_sites(self._sites[first:last+1]))

    def compute_WSSp(self):
        """
        Compute window summary statistics. Store data in ``_WSSp_stats``.
        """
        self._WSSp_stats = {}
        if len(self._params['WSSp']) == 0: return
        cs = egglib.stats.ComputeStats()
        cs.add_stats(*self._params['WSSp'])
        for k in self._struct_per_pop: self._WSSp_stats[k] = []
        for idx in self._selection:
            first = idx
            while first > 0:
                if self._chrom[first-1] != self._chrom[idx] or self._pos[first-1] < self._pos[idx]-self._params['wspan']:
                    break
                first -= 1
            last = idx
            while last < len(self._sites) - 1:
                if self._chrom[last+1] != self._chrom[idx] or self._pos[last+1] > self._pos[idx]+self._params['wspan']:
                    break
                last += 1
            for k, struct in self._struct_per_pop.iteritems():
                cs.configure(multi_hits=True, struct=struct)
                self._WSSp_stats[k].append(cs.process_sites(self._sites[first:last+1]))

    def export(self):
        """
        Export statistics value to output file.
        """
        try: f = open(self._params['output-file'], 'w')
        except IOError: raise ValueError, 'cannot open this file: {0}'.format(self._params['output-file'])
        f.write('ID')

        for ss in self._params['LSS']: f.write('\tLSS:{0}'.format(ss))
        for ss in self._params['LSSp']: f.write(''.join(['\tLSSp:{0}:{1}'.format(ss, p) for p in self._pops]))
        for ss in self._params['WSS']: f.write('\tWSS:{0}'.format(ss))
        for ss in self._params['WSSp']: f.write(''.join(['\tWSSp:{0}:{1}'.format(ss, p) for p in self._pops]))
        for ss in self._params['GSS']: f.write('\tGSS:{0}'.format(ss))
        for ss in self._params['GSSp']: f.write(''.join(['\tGSSp:{0}:{1}'.format(ss, p) for p in self._pops]))
        if self._params['SFS']:
            f.write('\t' + '\t'.join(['SFS:{0}'.format(i+1) for i in xrange(self._params['SFS-bins'])]))
        f.write('\n')

        GSS = '\t'.join(['NA' if self._GSS_stats[ss] is None else str(self._GSS_stats[ss]) for ss in self._params['GSS']])
        if len(GSS) > 0: GSS = '\t'+GSS
        GSSp = '\t'.join(['NA' if self._GSSp_stats[p][ss] is None else str(self._GSSp_stats[p][ss]) for ss in self._params['GSSp'] for p in self._pops])
        if len(GSSp) > 0: GSSp = '\t'+GSSp

        print '[summstats] exporting statistics for {0} loc{1}'.format(len(self._selection), 'us' if len(self._selection) == 1 else 'i')
        for idx, site_idx in enumerate(self._selection):
            f.write(self._id[site_idx])
            for ss in self._params['LSS']:
                if self._LSS_stats[idx][ss] is None: f.write('\tNA')
                else: f.write('\t{0}'.format(self._LSS_stats[idx][ss]))
            for ss in self._params['LSSp']:
                f.write('\t' + '\t'.join(['NA' if self._LSSp_stats[p][idx][ss] is None else str(self._LSSp_stats[p][idx][ss]) for p in self._pops]))
            for ss in self._params['WSS']:
                if self._WSS_stats[idx][ss] is None: f.write('\tNA')
                else: f.write('\t{0}'.format(self._WSS_stats[idx][ss]))
            for ss in self._params['WSSp']:
                f.write('\t' + '\t'.join(['NA' if self._WSSp_stats[p][idx][ss] is None else str(self._WSSp_stats[p][idx][ss]) for p in self._pops]))
            f.write(GSS)
            f.write(GSSp)
            if self._params['SFS']: f.write('\t' + '\t'.join(map(str, self._SFS)))
            f.write('\n')
        f.close()

if __name__ == '__main__':
    print '[summstats] starting'
    app = Application()

    try:
        app.run(sys.argv[1:])
    except Exception as e:
        sys.stderr.write('an error occurred: {0}\n'.format(e.message))
        print '[summstats] failed'
        if '_params' not in dir(app) or app._params['debug'] > 0:
            print '[summstat] full traceback'
            raise
    else:
        print '[summstats] done'



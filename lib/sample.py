#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a TCR repertoire sample
'''


class Sample():
    '''Represents a sample
    '''
    def __init__(self, name, clones=None, group=None):
        self.name = name
        self.clones = clones
        self.group = group

    def addclone(self, clone):
        self.clones.append(clone)

    def addclones(self, clones):
        self.clones.extend(clones)

    def setgroup(self, group):
        self.group = groupp




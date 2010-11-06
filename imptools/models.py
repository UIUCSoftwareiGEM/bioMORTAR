from django.db import models

class Organism(models.Model):
    ShortName=models.CharField(max_length=4, unique=True)
    FullName=models.CharField(max_length=100)   
    def __unicode__(self):    
        return self.ShortName
    class Meta:
        ordering = ['FullName']

class Compound(models.Model):
    KeggID=models.CharField(max_length=6)
    formula=models.CharField(max_length=20)
    Toxin=models.BooleanField()
    VAProduct=models.BooleanField()
    #Nontrivial=models.BooleanField()  
    def __unicode__(self):    
        return self.KeggID
    def get_names(self):
        return self.compoundname_set.all()
    def show_names(self):
        names=''
        for name in self.get_names():
            names=str(name) + '  ' + names
        return names
    show_names.short_description='Names'
    def MainName(self):
        return str(self.get_names()[0])
    class Meta:
        ordering = ['KeggID']
        
class CompoundName(models.Model):
    Name=models.CharField(max_length=40)
    comp=models.ForeignKey(Compound)
    def __unicode__(self):    
        return self.Name

    
    class Meta:    
        ordering = ['Name']

class Reactant(models.Model):
    reactname=models.ForeignKey(Compound, related_name='reactantof')
    stoich=models.CharField(max_length=2)
    def __unicode__(self):
        return self.reactname.KeggID
    class Meta:
                ordering = ['reactname']

class Product(models.Model):
    prodname=models.ForeignKey(Compound,related_name='productof')
    stoich=models.CharField(max_length=2)
    def __unicode__(self):
        return self.prodname.KeggID
    class Meta:
                ordering = ['prodname']

class Gene(models.Model):
    genename=models.CharField(max_length=30)
    organism=models.ForeignKey(Organism)
    def __unicode__(self):
        return self.genename
    class Meta:
        ordering = ['genename']

class Enzyme(models.Model):
    ECnumber=models.CharField(max_length=15)
    genes=models.ManyToManyField(Gene)
    def __unicode__(self):
        return self.ECnumber
    class Meta:
        ordering = ['ECnumber']
                
class Reaction(models.Model):
    KeggID=models.CharField(max_length=6)
    name=models.CharField(max_length=200)
    reversibility=models.BooleanField()
    reactants=models.ManyToManyField(Reactant)
    products=models.ManyToManyField(Product)
    organisms=models.ManyToManyField(Organism)
    enzymes=models.ManyToManyField(Enzyme)
    def __unicode__(self):
        return self.KeggID
    
    def atp(self):
        value=0
        #this returns the ATP use. Negative for consumpution, position for production
        for react in self.reactants.all():
            if str(react.reactname)=='C00002':
                value=value-int(react.stoich)
        for prod in self.products.all():
            if str(prod.prodname)=='C00002':
                value=value+int(react.stoich)
        return value
    def positive_atp(self,direction):
        value=0
        #this returns the ATP net consumption. positive value
        for react in self.reactants.all():
            if str(react.reactname)=='C00002':
                value=value-int(react.stoich)
        for prod in self.products.all():
            if str(prod.prodname)=='C00002':
                value=value+int(react.stoich)
        if direction>0:
            if value>=0:
                value=0
            else:
                value=-value
        else:
            if value<=0:
                value=0
            else:
                value=value   
        return value
    def enzymespresent(self):
        if len(self.enzymes.all() )>0:
            result=0
        else:
            result=1
        return result
    class Meta:
        ordering = ['KeggID']

import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_genes` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    # p[g][t]*p[t]=p[t][g]*p[g]
    # p[g][t]=p[t][g]*p[g]/p[t]

    jp=1
    no_gene=set(people.keys()) -one_gene -two_genes
    for person in no_gene:
        if person in have_trait:
            #determine pronb that person have trait and don't have gene
            #p(no ^ have)
            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                mom_pass_gene=(1-PROBS['mutation'])
                dad_pass_gene=(1-PROBS['mutation'])
                if mom in one_gene:
                    mom_pass_gene=PROBS['mutation']*(0.5)+0.5
                if dad in one_gene:
                    dad_pass_gene=PROBS['mutation']*(0.5)+0.5
                if mom in two_genes:
                    mom_pass_gene=PROBS['mutation'] 
                if dad in two_genes:
                    dad_pass_gene=PROBS['mutation']
                jp*=mom_pass_gene*dad_pass_gene*PROBS['trait'][0][True]
            else:
                jp*=PROBS['gene'][0]*PROBS['trait'][0][True]
        else:

            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                mom_pass_gene=(1-PROBS['mutation'])
                dad_pass_gene=(1-PROBS['mutation'])
                if mom in one_gene:
                    mom_pass_gene=PROBS['mutation']*(0.5)+0.5
                if dad in one_gene:
                    dad_pass_gene=PROBS['mutation']*(0.5)+0.5
                if mom in two_genes:
                    mom_pass_gene=PROBS['mutation'] 
                if dad in two_genes:
                    dad_pass_gene=PROBS['mutation']
                jp*=mom_pass_gene*dad_pass_gene*PROBS['trait'][0][False]
            else:
                jp*=PROBS['gene'][0]*PROBS['trait'][0][False]
    for person in one_gene:
        if person in have_trait:
            #determine pronb that person have trait and don't have gene
            #p(no ^ have)
            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                #mother gives gene and father does't give gene
                mul_f=0
                if mom in no_gene and dad in no_gene:
                    mul_f=PROBS['mutation']*(1-PROBS['mutation'])*2

                elif mom in no_gene and dad in one_gene:
                    mul_f=PROBS['mutation']*(0.5)+PROBS['mutation']*PROBS['mutation']*0.5+(1-PROBS['mutation'])*0.5*(1-PROBS['mutation'])

                elif mom in no_gene and dad in two_genes:
                    mul_f=PROBS['mutation']*PROBS['mutation']+(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in one_gene and dad in no_gene:
                    mul_f=(0.5)*(1-PROBS['mutation'])*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])

                elif mom in one_gene and dad in one_gene:
                    mul_f=((0.5)*(0.5)*(1-PROBS['mutation'])+(0.5)*(0.5)*(1-PROBS['mutation'])*(PROBS['mutation']))*2

                elif mom in one_gene and dad in two_genes:
                    mul_f=0.5*(1-PROBS['mutation'])*(PROBS['mutation'])+0.5*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in no_gene:
                    mul_f=(1-PROBS['mutation'])*(1-PROBS['mutation'])+PROBS['mutation']*PROBS['mutation']

                elif mom in two_genes and dad in one_gene:
                        mul_f=0.5*(1-PROBS['mutation'])*(PROBS['mutation'])+0.5*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in two_genes:
                    mul_f=(1-PROBS['mutation'])*(PROBS['mutation'])*2
                jp*=mul_f*PROBS['trait'][1][True]
            else:
                jp*=PROBS['gene'][1]*PROBS['trait'][1][True]

        else:
            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                #mother gives gene and father does't give gene
                mul_f=0
                if mom in no_gene and dad in no_gene:
                    mul_f=PROBS['mutation']*(1-PROBS['mutation'])*2

                elif mom in no_gene and dad in one_gene:
                    mul_f=PROBS['mutation']*(0.5)+PROBS['mutation']*PROBS['mutation']*0.5+(1-PROBS['mutation'])*0.5*(1-PROBS['mutation'])

                elif mom in no_gene and dad in two_genes:
                    mul_f=PROBS['mutation']*PROBS['mutation']+(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in one_gene and dad in no_gene:
                    mul_f=(0.5)*(1-PROBS['mutation'])*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])

                elif mom in one_gene and dad in one_gene:
                    mul_f=((0.5)*(0.5)*(1-PROBS['mutation'])+(0.5)*(0.5)*(1-PROBS['mutation'])*(PROBS['mutation']))*2

                elif mom in one_gene and dad in two_genes:
                    mul_f=0.5*(1-PROBS['mutation'])*(PROBS['mutation'])+0.5*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in no_gene:
                    mul_f=(1-PROBS['mutation'])*(1-PROBS['mutation'])+PROBS['mutation']*PROBS['mutation']

                elif mom in two_genes and dad in one_gene:
                        mul_f=0.5*(1-PROBS['mutation'])*(PROBS['mutation'])+0.5*(1-PROBS['mutation'])+0.5*(PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in two_genes:
                    mul_f=(1-PROBS['mutation'])*(PROBS['mutation'])*2
                jp*=mul_f*PROBS['trait'][1][False]
            else:
                jp*=PROBS['gene'][1]*PROBS['trait'][1][False]
    for person in two_genes:
        if person in have_trait:

            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                mul_f=0
                if mom in no_gene and dad in no_gene:
                    mul_f=(PROBS['mutation'])*(PROBS['mutation'])

                elif mom in no_gene and dad in one_gene:
                    mul_f=PROBS['mutation']*(0.5)*(1-PROBS['mutation'])
                elif mom in no_gene and dad in two_genes:
                    mul_f=PROBS['mutation']*(1-PROBS['mutation'])

                elif mom in one_gene and dad in no_gene:
                    mul_f=(0.5)*(1-PROBS['mutation'])*(PROBS['mutation'])

                elif mom in one_gene and dad in one_gene:
                    mul_f=(0.5)*(0.5)*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in one_gene and dad in two_genes:
                    mul_f=0.5*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in no_gene:
                    mul_f=(1-PROBS['mutation'])*(PROBS['mutation'])

                elif mom in two_genes and dad in one_gene:
                        mul_f=0.5*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in two_genes:
                    mul_f=(1-PROBS['mutation'])*(1-PROBS['mutation'])
                jp*=mul_f*(PROBS['trait'][2][True])
            else:
                jp*=PROBS['gene'][2]*PROBS['trait'][2][True]
        else:

            mom=people[person]['mother']
            dad=people[person]['father']
            if mom and dad:
                mul_f=0
                if mom in no_gene and dad in no_gene:
                    mul_f=(PROBS['mutation'])*(PROBS['mutation'])

                elif mom in no_gene and dad in one_gene:
                    mul_f=PROBS['mutation']*(0.5)*(1-PROBS['mutation'])
                elif mom in no_gene and dad in two_genes:
                    mul_f=PROBS['mutation']*(1-PROBS['mutation'])

                elif mom in one_gene and dad in no_gene:
                    mul_f=(0.5)*(1-PROBS['mutation'])*(PROBS['mutation'])

                elif mom in one_gene and dad in one_gene:
                    mul_f=(0.5)*(0.5)*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in one_gene and dad in two_genes:
                    mul_f=0.5*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in no_gene:
                    mul_f=(1-PROBS['mutation'])*(PROBS['mutation'])

                elif mom in two_genes and dad in one_gene:
                        mul_f=0.5*(1-PROBS['mutation'])*(1-PROBS['mutation'])

                elif mom in two_genes and dad in two_genes:
                    mul_f=(1-PROBS['mutation'])*(1-PROBS['mutation'])
                jp*=mul_f*(PROBS['trait'][2][False])
            else:
                jp*=PROBS['gene'][2]*PROBS['trait'][2][False]
    return jp 

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    no_gene=set(probabilities.keys())-one_gene-two_genes
    for person in one_gene:
        if person in have_trait:
            probabilities[person]['trait'][True]+=p
            probabilities[person]['gene'][1]+=p
        else:
            probabilities[person]['trait'][False]+=p
            probabilities[person]['gene'][1]+=p
    for person in two_genes:
        if person in have_trait:
            probabilities[person]['trait'][True]+=p
            probabilities[person]['gene'][2]+=p
        else:
            probabilities[person]['trait'][False]+=p
            probabilities[person]['gene'][2]+=p
    for person in no_gene:
        if person in have_trait:
            probabilities[person]['trait'][True]+=p
            probabilities[person]['gene'][0]+=p
        else:
            probabilities[person]['trait'][False]+=p
            probabilities[person]['gene'][0]+=p
    
    
def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        sm=probabilities[person]['gene'][0]+probabilities[person]['gene'][1]+probabilities[person]['gene'][2]
        probabilities[person]['gene'][0]/=sm 
        probabilities[person]['gene'][1]/=sm 
        probabilities[person]['gene'][2]/=sm 
        
        sm=probabilities[person]['trait'][True]+probabilities[person]['trait'][False]
        probabilities[person]['trait'][False]/=sm 
        probabilities[person]['trait'][True]/=sm
if __name__ == "__main__":
    main()

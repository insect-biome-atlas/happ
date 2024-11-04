# This code computes a distance between codons based
# on Grantham's distance (Grantham R, 1974: "Amino acid difference formula
# to help explain protein evolution")


aa_names <- c("Ser", "Arg", "Leu", "Pro", "Thr", "Ala", "Val", "Gly", "Ile", "Phe", "Tyr", "Cys", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Met", "Trp")

aa_dist <- matrix(nrow=20, ncol=20, byrow=TRUE, dimnames=list(aa_names,aa_names),
data=c(
0, 110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65, 80, 135, 177,
0, 0, 102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43, 86, 26, 96, 54, 91, 101, 
0, 0, 0, 98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113, 153, 107, 172, 138, 15, 61,
0, 0, 0, 0, 38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 91, 103, 108, 93, 87, 147,
0, 0, 0, 0, 0, 58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78, 85, 65, 81, 128,
0, 0, 0, 0, 0, 0, 64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107, 84, 148,
0, 0, 0, 0, 0, 0, 0, 109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88,
0, 0, 0, 0, 0, 0, 0, 0, 135, 153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184,
0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 33, 198, 94, 109, 149, 102, 168, 134, 10, 61,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22, 205, 100, 116, 158, 102, 177, 140, 28, 40,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 194, 83, 99, 143, 85, 160, 122, 36, 37,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 174, 154, 139, 202, 154, 170, 196, 215,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 68, 32, 81, 40, 87, 115,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46, 53, 61, 29, 101, 130,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 94, 23, 42, 142, 174,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 101, 56, 95, 110,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 160, 181,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126, 152,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
))

aa_dist <- aa_dist + t(aa_dist)
aa_dist <- aa_dist/(sum(aa_dist)/(20*20-20))    #normalize the distances


codon_names <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
                 "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
                 "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
                 "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
                 "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
                 "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
                 "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
                 "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")

# Translation table 5
aa_trans     <- c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser",
                 "Tyr", "Tyr", "Ter", "Ter", "Cys", "Cys", "Trp", "Trp",
                 "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro",
                 "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
                 "Ile", "Ile", "Met", "Met", "Thr", "Thr", "Thr", "Thr",
                 "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Ser", "Ser",
                 "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala",
                 "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly")

nuc_dist <- function(x,y) {
    dist <- 0
    for (i in 1:3) {
        c1 <- substring(x,i,i)
        c2 <- substring(y,i,i)

        if ((c1=="T" && c2=="C") || (c1=="C" && c2=="T"))
            dist <- dist + 1
        else if ((c1=="A" && c2=="G") || (c1=="G" && c2=="A"))
            dist <- dist + 1
        else if (c1!=c2)
            dist <- dist + 1
    }
    dist
}


codon_dist <- matrix(nrow=64,ncol=64,dimnames=list(codon_names,codon_names))
codon_aa_dist <- matrix(nrow=64,ncol=64,dimnames=list(codon_names,codon_names))
codon_nuc_dist <- matrix(nrow=64,ncol=64,dimnames=list(codon_names,codon_names))
codon_waa_dist <- matrix(nrow=64,ncol=64,dimnames=list(codon_names,codon_names))

for (i in 1:length(codon_names)) {
    codon1 <- codon_names[i]
    aa1 <- aa_trans[i]
    for (j in 1:length(codon_names)) {
        codon2 <- codon_names[j]
        aa2 <- aa_trans[j]

        codon_nuc_dist[i,j] <- nuc_dist(codon1, codon2)

        if (aa1=="Ter" || aa2=="Ter") {
          codon_dist[i,j] <- 1000
          codon_aa_dist[i,j] <- 1000
          codon_waa_dist[i,j] <- 1000
        } else if (aa1==aa2) {
          codon_dist[i,j] <- codon_nuc_dist[i,j]
          codon_aa_dist[i,j] <- 0
          codon_waa_dist[i,j] <- 0
        } else {
          codon_dist[i,j] <- codon_nuc_dist[i,j] + 10
          codon_aa_dist[i,j] <- 1
          codon_waa_dist[i,j] <- aa_dist[aa1,aa2]
        }
    }
}

# An arbitrary distance penalizing nonsynonymous changes
evo_dist <- function(x, y, start_pos=1) {

    x <- toupper(x)
    y <- toupper(y)
    
    dist <- 0

    for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {

        codon1 <- substring(x,i,i+2)
        codon2 <- substring(y,i,i+2)
        
        if (!grepl("-", codon1) && !grepl("-", codon2)) {
            dist <- dist + codon_dist[codon1,codon2]
        } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
                   (grepl("---", codon2) && !grepl("-", codon1))) {
          dist <- dist + 100
        }
    }

    dist / (nchar(x)-start_pos+1)
}

# A measure inspired by dN/dS metrics
dadn_dist <- function(x, y, start_pos=1) {
  
  x <- toupper(x)
  y <- toupper(y)
  
  da <- 0
  dn <- 0
  
  for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {
    
    codon1 <- substring(x,i,i+2)
    codon2 <- substring(y,i,i+2)

    if (!grepl("-", codon1) && !grepl("-", codon2)) {
      da <- da + codon_aa_dist[codon1,codon2]
      dn <- dn + codon_nuc_dist[codon1,codon2]
    } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
               (grepl("---", codon2) && !grepl("-", codon1))) {
      da <- da + 1
    }
      
  }
  
  da / (dn + 1E-10) # Guard against division by zero
}

# A measure inspired by dN/dS and also taking the biochemical
# simlarities of amino acids into account (a simple codon model
# of sorts)
wdadn_dist <- function(x, y, start_pos=1) {

  x <- toupper(x)
  y <- toupper(y)
  
  da <- 0
  dn <- 0
  
  for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {
    
    codon1 <- substring(x,i,i+2)
    codon2 <- substring(y,i,i+2)
    
    if (!grepl("-", codon1) && !grepl("-", codon2)) {
      da <- da + codon_waa_dist[codon1,codon2]
      dn <- dn + codon_nuc_dist[codon1,codon2]
    } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
               (grepl("---", codon2) && !grepl("-", codon1))) {
      da <- da + 1
    }
  }
  
  da / (dn + 1E-10) # Guard against division by zero
}


library(parallel)

if(!exists('kkp.kin')){
    kkp.kin<-readRDS('degree1_relatedness_king.rds')
    kkp.kin<-kkp.kin[kkp.kin$Kinship<0.354,]
}

if(!exists('age_sex')){
    age_sex<-readRDS('trio_info_age_sex.rds')
    age_sex$ids<-as.integer(rownames(age_sex))
}

find.trio<-function(){
    mat<-matrix(NA,nrow=nrow(age_sex),ncol=3)
    colnames(mat)<-c('Father','Mother','Child')
    mat[,3]<-age_sex$ids
    pmat.age<-matrix(NA,nrow=nrow(age_sex),ncol=2)
    colnames(pmat.age)<-c('Father','Mother')
    pmat.sex<-matrix(NA,nrow=nrow(age_sex),ncol=2)
    colnames(pmat.sex)<-c('Father','Mother')
    for(i in 1:nrow(kkp.kin)){
        x1<-kkp.kin[i,'ID1']
        x2<-kkp.kin[i,'ID2']
        age.1<-age_sex[age_sex$ids==x1,1]
        sex.1<-age_sex[age_sex$ids==x1,2]
        age.2<-age_sex[age_sex$ids==x2,1]
        sex.2<-age_sex[age_sex$ids==x2,2]
        if(age.1>=age.2){#x1 is parent #x2 is child
            ids<-mat[,3]==x2
            if(sum(ids)==1){
                if(all(is.na(mat[ids,1:2]))){#first meet id
                    if(sex.1==0){#1 is a mum
                        mat[ids,2]=x1
                        pmat.age[ids,2]=age.1
                        pmat.sex[ids,2]='F'
                    }
                    else{#1 is a dad
                        mat[ids,1]=x1
                        pmat.age[ids,1]=age.1
                        pmat.sex[ids,1]='M'
                    }
                }
                else if(all(!is.na(mat[ids,1:2]))){#both been met,add a row, an odd case
                    if(sex.1==0){#1 is a mum
                        m.tmp<-c(NA,x1,x2)
                        age.tmp<-c(NA,age.1)
                        sex.tmp<-c(NA,'F')
                    }
                    else{#1 is a dad
                        m.tmp<-c(x1,NA,x2)
                        age.tmp<-c(age.1,NA)
                        sex.tmp<-c('M',NA)
                    }
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
                else if(!is.na(mat[ids,1]) && is.na(mat[ids,2])){# previously found to be duo: Father and child
                    if(sex.1==0){#1 is a mum
                        mat[ids,2]=x1
                        pmat.age[ids,2]=age.1
                        pmat.sex[ids,2]='F'
                    }
                    else{#two mums?
                        m.tmp<-c(x1,NA,x2)
                        age.tmp<-c(NA,age.1)
                        sex.tmp<-c(NA,'M')
                        mat<-rbind(mat,m.tmp)
                        pmat.age<-rbind(pmat.age,age.tmp)
                        pmat.sex<-rbind(pmat.sex,sex.tmp)
                    }
                }
                else{#duo, mum and child
                    if(sex.1==1){#1 is a dad
                        mat[ids,1]=x1
                        pmat.age[ids,1]=age.1
                        pmat.sex[ids,1]='M'
                    }
                    else{#1 is another mum
                        m.tmp<-c(NA,x1,x2)
                        age.tmp<-c(NA,age.1)
                        sex.tmp<-c(NA,'F')
                        mat<-rbind(mat,m.tmp)
                        pmat.age<-rbind(pmat.age,age.tmp)
                        pmat.sex<-rbind(pmat.sex,sex.tmp)
                    }
                }
            }
            else{#odd things happen,just append
                if(sex.1==1){#1 is a dad
                    m.tmp<-c(x1,NA,x2)
                    age.tmp<-c(age.1,NA)
                    sex.tmp<-c('F',NA)
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
                else{#1 is a mum
                    m.tmp<-c(NA,x1,x2)
                    age.tmp<-c(NA,age.1)
                    sex.tmp<-c(NA,'F')
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
            }
        }
        else{
            ids<-mat[,3]==x1
            if(sum(ids)==1){
                if(all(is.na(mat[ids,1:2]))){#first meet id
                    if(sex.2==0){#2 is a mum
                        mat[ids,2]=x2
                        pmat.age[ids,2]=age.2
                        pmat.sex[ids,2]='F'
                    }
                    else{#2 is a dad
                        mat[ids,1]=x2
                        pmat.age[ids,1]=age.2
                        pmat.sex[ids,1]='M'
                    }
                }
                else if(all(!is.na(mat[ids,1:2]))){#both been met,add a row, an odd case
                    if(sex.2==0){#2 is a mum
                        m.tmp<-c(NA,x2,x1)
                        age.tmp<-c(NA,age.2)
                        sex.tmp<-c(NA,'F')
                    }
                    else{#2 is a dad
                        m.tmp<-c(x2,NA,x1)
                        age.tmp<-c(age.2,NA)
                        sex.tmp<-c('M',NA)
                    }
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
                else if(!is.na(mat[ids,1]) && is.na(mat[ids,2])){# previously found to be duo: Father and child
                    if(sex.2==0){#2 is a mum
                        mat[ids,2]=x2
                        pmat.age[ids,2]=age.2
                        pmat.sex[ids,2]='F'
                    }
                    else{#two mums?
                        m.tmp<-c(x2,NA,x1)
                        age.tmp<-c(NA,age.2)
                        sex.tmp<-c(NA,'M')
                        mat<-rbind(mat,m.tmp)
                        pmat.age<-rbind(pmat.age,age.tmp)
                        pmat.sex<-rbind(pmat.sex,sex.tmp)
                    }
                }
                else{#duo, mum and child
                    if(sex.2==1){#1 is a dad
                        mat[ids,1]=x2
                        pmat.age[ids,1]=age.2
                        pmat.sex[ids,1]='M'
                    }
                    else{#2 is another mum
                        m.tmp<-c(NA,x2,x1)
                        age.tmp<-c(NA,age.2)
                        sex.tmp<-c(NA,'F')
                        mat<-rbind(mat,m.tmp)
                        pmat.age<-rbind(pmat.age,age.tmp)
                        pmat.sex<-rbind(pmat.sex,sex.tmp)
                    }
                }
            }
            else{#odd things happen,just append
                if(sex.2==1){#2 is a dad
                    m.tmp<-c(x2,NA,x1)
                    age.tmp<-c(age.2,NA)
                    sex.tmp<-c('F',NA)
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
                else{#1 is a mum
                    m.tmp<-c(NA,x2,x1)
                    age.tmp<-c(NA,age.2)
                    sex.tmp<-c(NA,'F')
                    mat<-rbind(mat,m.tmp)
                    pmat.age<-rbind(pmat.age,age.tmp)
                    pmat.sex<-rbind(pmat.sex,sex.tmp)
                }
            }
        
        }
    }
    return(list(mat,pmat.age,pmat.sex))
}

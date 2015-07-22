get.kmer.events <- function(f5) {

        ev <- get.events(f5)

        tmp.events <- list(template  = ev$template,
                           complement= ev$complement)

        readAln <- get.align(f5)

        readAln$template <- as.numeric(readAln$template)+1;
        readAln$complement <- as.numeric(readAln$complement)+1;
        readAln$template[readAln$template == 0] <- NA;
        readAln$complement[readAln$complement == 0] <- NA;

         for(type in c("template","complement")){
                for(stat in c("mean","stdv")){
                        readAln[[paste(stat,type,sep=".")]] <-
                                as.numeric((tmp.events[[type]])[[stat]][readAln[[type]]]);
                }
                for(stat in c("start","length")){
                        readAln[[paste(stat,type,sep=".")]] <-
                                as.numeric((tmp.events[[type]])[[stat]][readAln[[type]]]) * 1;
                }

                for(stat in c("model_state","p_model_state","mp_state","p_mp_state")){
                        readAln[[paste(stat,type,sep=".")]] <-
                                as.vector((tmp.events[[type]])[[stat]][readAln[[type]]]);
                }
        }
        #readAln$move <- kmer2move(readAln$kmer);
        #readAln$called.pos <- cumsum(readAln$move) + 1;
        rownames(readAln) <- 0:(dim(readAln)[1]-1);

        return(readAln)

}



get.align <- function(f5) {

        fid <- H5Fopen(f5)

        gid <- H5Gopen(fid, "/Analyses/Basecall_2D_000/BaseCalled_2D/");

        did <- H5Dopen(gid, "Alignment")

        aln <- (H5Dread(did, bit64conversion = "bit64"))

        H5Dclose(did)
        H5Gclose(gid)
        H5Fclose(fid)

        return(aln)
}



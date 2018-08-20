#' Build index for mapping using Salmon
#'
#' \code{build_index} for mapping reads using Salmon.
#'
#' @inheritParams run_dea
#' @param kmer k-mer size for indexing. default is 31. See \code{'Salmon'}
#' for details.
#' @param ens_release version of Ensembl release.
#' @param dest.dir directory where all the files will be saved.
#'
#' @return directory of index files
#'
#' @references
#' Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and
#' Carl Kingsford (2017): Salmon provides fast and bias-aware
#' quantification of transcript expression. Nature methods, 14(4),417.
#' \url{https://www.nature.com/articles/nmeth.4197}
#'
#' @export

build_index <- function(species = c("human", "mouse"), kmer = 31,
                        ens_release = 93, dest.dir) {

  species <- match.arg(species, c("human", "mouse"))
  create_dir(dest.dir)

  if (species == "human") {
    run_cmd(paste0("wget -O ", dest.dir, "/Homo_sapiens.GRCh38.cdna.",
                  "all.fa.gz ftp://ftp.ensembl.org/pub/release-",
                  ens_release, "/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))
    run_cmd(paste0("wget -O ", dest.dir,"/Homo_sapiens.GRCh38.ncrna.fa.gz",
                  " ftp://ftp.ensembl.org/pub/release-", ens_release,
                  "/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"))
    run_cmd(paste0("gunzip -c ", dest.dir, "/Homo_sapiens.GRCh38.cdna.all.",
                  "fa.gz ", dest.dir, "/Homo_sapiens.GRCh38.ncrna.fa.gz > ", dest.dir,
                  "/Homo_sapiens.GRCh38.release", ens_release, ".cdna.ncrna.fa"))
    run_cmd(paste0("salmon index -t ", dest.dir, "/Homo_sapiens.GRCh38.",
                  "release", ens_release, ".cdna.ncrna.fa -i ", dest.dir,
                  "/human_transcripts_release", ens_release, "_index"))
  }
  if (species == "mouse") {
    run_cmd(paste0("wget -O ", dest.dir, "/Mus_musculus.GRCm38.cdna.",
                  "all.fa.gz ftp://ftp.ensembl.org/pub/release-",
                  ens_release, "/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"))
    run_cmd(paste0("wget -O ", dest.dir, "/Mus_musculus.GRCm38.ncrna.fa.gz",
                  " ftp://ftp.ensembl.org/pub/release-", ens_release,
                  "/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"))
    run_cmd(paste0("gunzip -c ", dest.dir, "/Mus_musculus.GRCm38.cdna.all.",
                  "fa.gz ", dest.dir, "/Mus_musculus.GRCm38.ncrna.fa.gz > ", dest.dir,
                  "/Mus_musculus.GRCm38.release", ens_release, ".cdna.ncrna.fa"))
    run_cmd(paste0("salmon index -t ", dest.dir, "/Mus_musculus.GRCm38.",
                  "release", ens_release, ".cdna.ncrna.fa -i ", dest.dir,
                  "/mouse_transcripts_release", ens_release,"_index"))
  }
}



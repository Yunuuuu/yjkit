curl_downloads <- function(urls, destfiles, mode) {
    curl_pool <- curl::new_pool()
    .mapply(
        function(url, destfile) {
            force(url)
            force(destfile)
            curl::curl_fetch_multi(
                url = url, done = function(req) {
                    switch(mode,
                        w = data.table::fwrite(
                            list(rawToChar(req$content)), 
                            destfile, sep = ""
                        ),
                        wb = writeBin(req$content, destfile)
                    )
                }, fail = function(msg) {
                    rlang::abort(c(paste0("Error: ", url), msg))
                },
                pool = curl_pool,
                handle = curl::new_handle(
                    timeout = 120L, connecttimeout = 60L,
                    failonerror = TRUE
                )
            )
        },
        list(urls, destfiles),
        MoreArgs = NULL
    )
    curl::multi_run(pool = curl_pool)
    invisible(destfiles)
}

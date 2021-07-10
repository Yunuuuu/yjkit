#' @importFrom ggplot2 aes
#' @export
ggplot2::aes

#' Sequential, diverging and qualitative colour scales from \code{paletteer}
#'
#' @description  The choices of color palettes in R can be quite overwhelming
#'   with palettes spread over many packages with many different API's. This
#'   packages aims to collect all color palettes across the R ecosystem under
#'   the same package with a streamlined API. See
#'   \code{\link[paletteer]{paletteer-package}}
#' @param palette Name of palette as a string. Can be in the form of
#'   \code{packagename::palettename}. Details see
#'   \href{https://emilhvitfeldt.github.io/paletteer/}{\code{paletteer}}. for
#'   \code{palette} in \code{ggsci}, we will use regular expression to match. If
#'   match nothing in \code{paletteer palette}.
#'   \code{\link[ggplot2:scale_manual]{scale_discrete_manual}} will be used.
#' @param scale the aesthetic mapping the palette
#' @param type the type of scale, a scalar character, can be "\code{discrete}",
#'   "\code{continuous}", and "\code{binned}". Default: "\code{discrete}"
#' @param direction Sets the order of colours in the scale. If 1, the default,
#'   colours are as output by
#'   \code{\link[paletteer:ggplot2-scales-continuous]{paletteer}}. If -1, the
#'   order of colours is reversed. Default: \code{1}.
#' @param ... Additional arguments pass on to
#'   \code{\link[ggplot2:discrete_scale]{discrete scale}},
#'   \code{\link[ggplot2:scale_gradient]{gradient continuous scale}} or
#'   \code{\link[ggplot2:scale_steps]{stepsn scale}}.
#' @return  A Scale object that can be added to a ggplot object
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   ggscale_paletteer()
#' @export
ggscale_paletteer <- function(palette = "nejm",
                              scale = "colour",
                              type = c("discrete", "continuous", "binned"),
                              direction = 1,
                              ...){

  if (!rlang::is_scalar_character(scale)){

    warning(paste0("scale should be a scalar character, ",
                   "we'll take the first one and as.character"))
    scale <- as.character(scale)[[1]]

  }

  if (!direction %in% c(-1, 1)) stop(
    "direction must in c(-1, 1)"
  )

  scale_match_color <- pmatch(scale, "color", nomatch = 0)

  if (scale_match_color == 1){
    scale <- "color"
  }

  if (scale_match_color == 0){
    scale <- match.arg(scale, c("colour", "fill"))
  }

  if (rlang::is_scalar_character(palette)) {

    type <- match.arg(type)

    type <- switch(
      type,
      discrete = "d",
      continuous = "c",
      binned = "binned"
    )

    palette_fn <- rlang::eval_tidy( rlang::call2(
      "::", rlang::expr(paletteer), rlang::sym(
        stringr::str_c("scale_", scale, "_paletteer_", type)
      )
    ) )

    all_paletteer_palette <- stringr::str_c(
      c(paletteer::palettes_d_names$package,
        paletteer::palettes_c_names$package,
        paletteer::palettes_dynamic_names$palette),
      c(paletteer::palettes_d_names$palette,
        paletteer::palettes_c_names$palette,
        paletteer::palettes_dynamic_names$palette),
      sep = "::"
    )

    ggsci_palette <- stringr::str_subset(
      stringr::str_subset(all_paletteer_palette, "^ggsci::"),
      pattern = stringr::str_c("ggsci::.*", palette, sep = "")
    )

    if (length(ggsci_palette) > 0) {

      if (length(ggsci_palette) > 1) {
        warning("more than one ggsci palette found, we'll take the first one")
        palette <- ggsci_palette[[1]]
      } else {
        palette <- ggsci_palette
      }

      return(palette_fn(palette = palette,
                        direction = direction,
                        ...))

    } else if ( palette %in% all_paletteer_palette ){

      return(palette_fn(palette = palette,
                        direction = direction,
                        ...))

    }
  }

  if (direction == -1) palette <- rev(palette)
  return( ggplot2::scale_discrete_manual(
    aesthetics = scale,
    values = palette, ...
  ) )

}

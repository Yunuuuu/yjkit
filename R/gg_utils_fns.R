#' @importFrom ggplot2 aes
#' @export
ggplot2::aes

# valid_text_aes <- function(x, type){
#
#   if (is.null(names(x))) stop(x, " must has names")
#
#   element_text_to_gpar <- c(
#     family = "fontfamily",
#     face = "fontface",
#     color = "col",
#     colour = "col",
#     size = "fontsize",
#     angle = "rot",
#     lineheight = "lineheight",
#     margin = "margin"
#   )
# }

#' Sequential, diverging and qualitative colour scales from \code{paletteer}
#'
#' @description  The choices of color palettes in R can be quite overwhelming
#'   with palettes spread over many packages with many different API's. This
#'   packages aims to collect all color palettes across the R ecosystem under
#'   the same package with a streamlined API. See
#'   \code{\link[paletteer]{paletteer-package}}
#' @param palette Name of palette as a string. Must be in the form of
#'   \code{packagename::palettename}. Details see
#'   \code{\href{https://emilhvitfeldt.github.io/paletteer/}{paletteer}}. for
#'   \code{palette} in \code{ggsci}, we will use regular expression to match.
#' @param scale the aesthetic mapping the palette
#' @param type the type of scale, a scalar character, can be "\code{discrete}",
#'   "\code{continuous}", and "\code{binned}". Default: "\code{discrete}"
#' @param direction Sets the order of colours in the scale. If 1, the default,
#'   colours are as output by
#'   \code{\link[paletteer:scale_color_paletteer_c]{paletteer}}. If -1, the
#'   order of colours is reversed. Default: \code{1}.
#' @param ... Additional arguments pass on to
#'   \code{\link[ggplot2:discrete_scale]{discrete scale}},
#'   \code{\link[ggplot2:scale_colour_gradientn]{gradient continuous scale}} or
#'   \code{\link[ggplot2:scale_colour_stepsn]{stepsn scale}}.
#' @return  A Scale object that can be added to a ggplot object
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   gg_scale_paletteer()
ggscale_paletteer <- function(palette = "nejm",
                              scale = c("colour", "color", "fill"),
                              type = c("discrete", "continuous", "binned"),
                              direction = 1,
                              ...){

  scale <- match.arg(scale)
  type <- match.arg(type)

  type <- switch(
    type,
    discrete = "d",
    continuous = "c",
    binned = "binned"
  )

  palette_fn <- rlang::eval_tidy( rlang::expr(
    `::`(paletteer, !!rlang::sym(
      stringr::str_c("scale_", scale, "_paletteer_", type)
    ))
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
    pattern = stringr::str_c("ggsci::.*", palette)
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

  } else {

    return( ggplot2::scale_color_manual( values = palette, ... ) )

  }

}


# annotate_npc <- function(label, x, y, ...)
# {
#   ggplot2::annotation_custom(grid::textGrob(
#     x = grid::unit(x, "npc"),
#     y = grid::unit(y, "npc"), label = label, ...))
# }

#' Create an annotation layer
#'
#' This function adds geoms to a plot with x and y mapping to Normalised Parent
#' Coordinates
#'
#' @param geom character Name of geom to use for annotation. Can be one of
#'   \code{"grob"}, \code{"plot"}, \code{"text"}, \code{"label"}, and
#'   \code{"table"}. Details see \code{\link[ggpmisc]{annotate}}
#' @param x,y the position of added label in Normalised Parent
#'   Coordinates(\code{npc})
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}.
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @details See \code{\link{ggplot2::annotate}}
#' @export
ggannotate_npc <- function(geom, x, y, na.rm = FALSE, ...){

  if (!geom %in% c("text", "label")) stop(
    'geom must be one of "text" and "label"'
  )
  ggplot2::annotate(geom = stringr::str_c(geom, "_npc", sep = ""),
                    x_npc = x,
                    y_npc = y,
                    na.rm = na.rm,
                    ...)
}

#' Create an ggplot2 object with annotation
#'
#' This function Creates an ggplot2 object with annotation.
#'
#' @param data Default dataset to use for plot. If not specified, must be
#'   supplied in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot. If not
#'   specified, must be supplied in each layer added to the plot.
#' @param label the character which will be added to the plot, will be connected
#'   by \code{stringr::str_c(names(label), label), collapse = label_sep, sep =
#'   ": "}
#' @param format_label if \code{TRUE}, will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length 2 or a character vector \code{(one of "title",
#'   "subtitle", "caption", and "tag")} of length one.
#' @param label_justification label_justification must be a element-two (or one)
#'   character or numeric vector. Possible string values are: "left", "right",
#'   "centre", "center", "bottom", and "top". For numeric values, 0 means left
#'   (bottom) alignment and 1 means right (top) alignment. Details see
#'   \code{\link[grid]{valid.just}}. Default: for numeric \code{label_position},
#'   default \code{label_justification} is \code{c(0, 1.15)}, for character
#'   \code{label_position}, default \code{label_justification} is \code{c(0, 0)}
#' @param label_sep see argument \code{label}. String to insert between each of
#'   the \code{label}.
#' @param ggtheme a ggplot2 theme object. Provided as a ggplot2 theme object or
#'   a function (can be a function name as a string) implemented as a ggplot2
#'   theme object.
#' @param ... additional arguments control the \code{label}. See
#'   \code{\link[grid]{textGrob}} and \code{\link[ggplot2]{element_text}}
#' @return a ggplot2 object with annotation
#' @examples
#'   gganno(mtcars, aes(cyl, mpg), label = "a")
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
gganno_text <- function(data, mapping = NULL,
                        label, format_label = TRUE,
                        label_position = c(0.02, 1),
                        label_justification = NULL,
                        label_sep = NULL, ggtheme = ggthemes::theme_few(),
                        ...){

  if (!class(label_position) %in% c("numeric", "character")) {
    stop(
      paste0("label_position must be a element-two numeric ",
             "or a element-one character (",
             'one of "title", "subtitle", "caption", and "tag" (length 1)',
             ") vector")
    )
  }

  if (is.null(label_justification)) label_justification <- switch(
    class(label_position),
    numeric = c(0, 1.15),
    character = c(0, 0)
  )

  if (!class(label_justification) %in% c("numeric", "character")) {
    stop(paste0(
      "label_justification must be a element-two (or one) character ",
      'or numeric vector. Possible string values are: "left", ',
      '"right", "centre", "center", "bottom", and "top". For numeric values, ',
      "0 means left (bottom) alignment and 1 means right (top) alignment."
    ))
  }

  label_justification <- grid::valid.just(label_justification)

  if (!ggplot2::is.theme(ggtheme)) {
    ggtheme <- rlang::exec(ggtheme)
  }

  if (!ggplot2::is.theme(ggtheme)) stop(
    "ggtheme must be a ggplot2 theme object or a function (can be a function name as a string) implemented as a ggplot2 theme object"
  )
  if (is.null(mapping)) mapping <- aes()

  res <- ggplot2::ggplot(data = data, mapping = mapping) + ggtheme

  if (format_label & is.numeric(label)) {
    label_format <- format_num(label)
  } else {
    label_format <- as.character(label)
  }

  if (is.numeric(label_position)) {

    if (length(label_position) != 2) stop("the numeric vector label_position must has a length of 2")

    if (is.null(label_sep)) label_sep = "\n"

    res <- res + ggannotate_npc(
      geom = "text",
      label = stringr::str_c(names(label), label_format, sep = ": ",
                             collapse = label_sep),
      x = label_position[[1]],
      y = label_position[[2]],
      hjust = label_justification[[1]],
      vjust = label_justification[[2]],
      ...
    )

  }

  if (is.character(label_position)) {

    label_position <- label_position[[1]]

    if ( label_position %in% c("title", "subtitle", "caption", "tag") ) stop(
      paste0('label_position must in c("title", "subtitle", "caption", "tag")')
    )

    if (is.null(label_sep)) label_sep <- "; "

    res + ggplot2::labs(
      !!label_position := stringr::str_c(names(label),
                                         label_format, sep = ": ",
                                         collapse = label_sep)
    ) + ggplot2::theme(
      !!stringr::str_c("plot.", label_position) := ggplot2::element_text(
        hjust = label_justification[[1]],
        vjust = label_justification[[2]],
        ...
      )
    )

  }

  res

}



#' Compute npc coordinates
#'
#' Convert character-encoded positions to npc units and keep numeric as npc
#' values.
#'
#' @param x numeric or character vector of coordinates.
#' @return A numeric vector with values representing npc coordinates.
valid_npc <- function(x){

  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) return(
    unname(c(left = 0, center = 0.5, right = 1,
             bottom = 0, middle = 0.5, top = 1)[x]))

  if (is.numeric(x)) {
    return(x)
  } else {
    return(as.numeric(x))
  }

}

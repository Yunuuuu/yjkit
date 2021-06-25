#' @export
#' @rdname geom_text_npc
#' @param label.padding Amount of padding around label. Defaults to 0.25 lines.
#' @param label.r Radius of rounded corners. Defaults to 0.15 lines.
#' @param label.size Size of label border, in mm.
geom_label_npc <- function(mapping = NULL, data = NULL,
                           stat = "identity", position = "identity",
                           ...,
                           parse = FALSE,
                           nudge_x = 0,
                           nudge_y = 0,
                           label.padding = grid::unit(0.25, "lines"),
                           label.r = grid::unit(0.15, "lines"),
                           label.size = 0.25,
                           na.rm = FALSE,
                           show.legend = NA,
                           inherit.aes = TRUE) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      abort("You must specify either `position` or `nudge_x`/`nudge_y`.")
    }

    position <- ggplot2::position_nudge(nudge_x, nudge_y)
  }

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomLabelNpc,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      parse = parse,
      label.padding = label.padding,
      label.r = label.r,
      label.size = label.size,
      na.rm = na.rm,
      ...
    )
  )
}

#' Base ggproto classes
#'
#' This geom is identical to 'ggplot2' \code{\link[ggplot2]{geom_text}}
#'   and  except that it interprets \code{x} and \code{y} positions in
#'   \code{npc} units. It translates \code{x} and \code{y} coordinates from npc
#'   units to native data units.
#'
#' @details
#'   See \code{\link[ggplot2]{ggplot2-ggproto}}
#' @seealso
#'   \code{\link[ggplot2]{ggproto}}
#' @rdname ggplot2-ggproto-npc
#' @details
#' @format NULL
#' @usage NULL
#' @export
GeomLabelNpc <- ggplot2::ggproto(
  "GeomLabelNpc", ggplot2::Geom,
  required_aes = c("x_npc", "y_npc", "label"),

  default_aes = ggplot2::GeomLabel$default_aes,

  draw_panel = function(data, panel_params, coord, parse = FALSE,
                        na.rm = FALSE,
                        label.padding = grid::unit(0.25, "lines"),
                        label.r = grid::unit(0.15, "lines"),
                        label.size = 0.25) {

    data$x_npc <- grid::valid.just(data$x_npc)
    data$y_npc <- grid::valid.just(data$y_npc)

    ranges <- coord$backtransform_range(panel_params)

    data$x <- ranges$x[1] + data$x_npc * (ranges$x[2] - ranges$x[1])
    data$y <- ranges$y[1] + data$y_npc * (ranges$y[2] - ranges$y[1])

    ggplot2::GeomLabel$draw_panel(data = data,
                                  panel_params = panel_params,
                                  coord = coord,
                                  parse = parse,
                                  na.rm = na.rm,
                                  label.padding = label.padding,
                                  label.r = label.r,
                                  label.size = label.size)
  },

  draw_key = ggplot2::GeomLabel$draw_key
)

#' Create an annotation layer with text or label annotation
#'
#' This function Creates an ggplot2 object with annotation.
#'
#' @param geom character Name of geom to use for annotation. Can be one of
#'   \code{"text"} and \code{"label"}.
#' @param label the character which will be added to the plot, will be connected
#'   by \code{stringr::str_c(names(label), label), collapse = label_sep, sep =
#'   ": "}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length two or a character vector with length one (one of
#'   \code{"lefttop", "righttop", "leftbottom", and "rightbottom"}) or length
#'   two (See below \code{label_justification}). Default: \code{"lefttop"}.
#' @param label_justification label_justification must be a element-two (or one)
#'   character or numeric vector. Possible string values are: "left", "right",
#'   "centre", "center", "bottom", and "top". For numeric values, 0 means left
#'   (bottom) alignment and 1 means right (top) alignment.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These
#'   are often aesthetics, used to set an aesthetic to a fixed value, like
#'   colour = "red" or size = 3. They may also be parameters to the paired
#'   geom/stat. See \code{\link[ggplot2:geom_text]{geom_text}} and
#'   \code{\link[ggplot2:geom_text]{geom_label}}.
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @return a ggplot2 layer object with annotation
#' @examples
#' ggannotate_npc(label = "a", label_position = c(0.5, 0.5))
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
ggannotate_npc <- function(geom = "text",
                           label, label_position = "lefttop",
                           label_justification = NULL,
                           ..., na.rm = FALSE) {
  if (is.numeric(label_position)) {
    if (!identical(length(label_position), 2L)) {
      stop("label_position should be a length-two when it is ",
        "a numeric vector",
        call. = FALSE
      )
    }

    if (is.null(label_justification)) {
      label_justification <- c(0.5, 0.5)
    }
  } else if (is.character(label_position)) {
    if (identical(length(label_position), 1L)) {
      if (label_position %in% c(
        "lefttop", "righttop", "leftbottom", "rightbottom"
      )) {
        if (is.null(label_justification)) {
          label_justification <- switch(label_position,
            lefttop = c(0, 0),
            righttop = c(1, 0),
            leftbottom = c(0, 1),
            rightbottom = c(1, 1)
          )
        }

        label_position <- switch(label_position,
          lefttop = c(0.02, 1.02),
          righttop = c(0.98, 1.02),
          leftbottom = c(0.02, -0.02),
          rightbottom = c(0.98, -0.02)
        )
      } else {
        stop("label_position with a length-one character vector should be ",
          "one of ",
          '"lefttop", "righttop", "leftbottom", and "rightbottom".',
          call. = FALSE
        )
      }
    } else if (identical(length(label_position), 2L)) {
      if (!label_position[[1]] %in% c("left", "right", "center", "middle")) {
        stop("the first element of label_position should be in ",
          "one of ", 'c("left", "right", "center", "middle").',
          call. = FALSE
        )
      }

      if (!label_position[[2]] %in% c("bottom", "top", "center", "middle")) {
        stop("the second element of label_position should be in ",
          "one of ", 'c("bottom", "top", "center", "middle").',
          call. = FALSE
        )
      }
    } else {
      stop("label_position should be a ",
        "length-one or length-two when it is a character vector",
        call. = FALSE
      )
    }

    label_position <- valid_npc(label_position)

    if (is.null(label_justification)) {
      label_justification <- c(0.5, 0.5)
    }
  } else {
    stop("label_position should be a length-two numeric vector ",
      "or a length-one or -two character vector.",
      call. = FALSE
    )
  }

  ggannotate_npc_helper(
    geom = geom,
    label = label,
    x = label_position[[1]],
    y = label_position[[2]],
    hjust = label_justification[[1]],
    vjust = label_justification[[2]],
    ...,
    na.rm = na.rm
  )
}

# annotate_npc <- function(label, x, y, ...)
# {
#   ggplot2::annotation_custom(grid::textGrob(
#     x = grid::unit(x, "npc"),
#     y = grid::unit(y, "npc"), label = label, ...))
# }

#' Create an ggplot2 annotation layer
#'
#' This function adds geoms to a plot with x and y mapping to Normalised Parent
#' Coordinates
#'
#' @param geom character Name of geom to use for annotation. Can be one of
#'   \code{"text"} and \code{"label"}.
#' @param x,y the position of added label in Normalised Parent
#'   Coordinates(\code{npc})
#' @param na.rm If \code{FALSE} (the default), removes missing values with a
#'   warning.  If \code{TRUE} silently removes missing values.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These
#'   are often aesthetics, used to set an aesthetic to a fixed value, like
#'   colour = "red" or size = 3. They may also be parameters to the paired
#'   geom/stat. See \code{\link[ggplot2:geom_text]{geom_text}} and
#'   \code{\link[ggplot2:geom_text]{geom_label}}.
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @details See \code{\link[ggplot2]{annotate}}
ggannotate_npc_helper <- function(geom, x = NULL, y = NULL,
                                  ..., na.rm = FALSE) {
  if (!geom %in% c("text", "label")) {
    stop(
      'geom must be one of "text" and "label"',
      call. = FALSE
    )
  }
  position <- purrr::compact(list(
    x_npc = x, y_npc = y
  ))
  aesthetics <- c(position, list(...))

  # Check that all aesthetic have compatible lengths
  lengths <- vapply(aesthetics, length, integer(1))
  n <- unique(lengths)

  # if there is more than one unique length, ignore constants
  if (length(n) > 1L) {
    n <- setdiff(n, 1L)
  }

  # if there is still more than one unique length, we error out
  if (length(n) > 1L) {
    bad <- lengths != 1L
    details <- paste(names(aesthetics)[bad], " (", lengths[bad], ")",
      sep = "", collapse = ", "
    )
    rlang::abort(glue::glue("Unequal parameter lengths: {details}"))
  }

  data <- vctrs::new_data_frame(position, n = n)
  ggplot2::layer(
    geom = stringr::str_c(geom, "_npc", sep = ""),
    params = list(
      na.rm = na.rm,
      ...
    ),
    stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity,
    data = data,
    mapping = ggplot2::aes_all(names(data)),
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}

#' Compute npc coordinates
#'
#' Convert character-encoded positions to npc units and keep numeric as npc
#' values.
#'
#' @param x numeric or character vector of coordinates.
#' @return A numeric vector with values representing npc coordinates.
valid_npc <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    return(
      unname(c(
        left = 0, center = 0.5, right = 1,
        bottom = 0, middle = 0.5, top = 1
      )[x])
    )
  }

  if (is.numeric(x)) {
    return(x)
  } else {
    return(as.numeric(x))
  }
}

setMethod("cleave", "Proteins",
          function(x, ...) {
              x@pfeatures <- cleave(aa(x), ...)
              x
          })

setGeneric('SubsetPoints',
           function(BITraster, points, variables, years, parallel = FALSE) standardGeneric('SubsetPoints'))
setGeneric('TemplateBIT',
           function(template, dir, variables, years) standardGeneric('TemplateBIT'))
setGeneric('MakeMask',
           function(BITraster, mask, filename) standardGeneric('MakeMask'))

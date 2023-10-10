# PointedSDMs 1.3.1

#### Changes and fixes since previous version:

-   Changed *SetophagaData* to an s*f* object.
-   Fixed issue regarding *crs* object coming from *sf.*
-   Fixed the multinomial model for when count data is considered.
-   Changed the speciesSpatial model: instead of being logical, the variable may now take on the values *NULL* (for no species specific spatial effects), *copy* to create a copy model of the spatial effect for the species across the dataset, or *individual* for creating independent spatial effects for each species.
-   Added additional checks to fix the errors regarding creating a *sf* polygon object from the mesh.
-   Documentation and spelling fixes.

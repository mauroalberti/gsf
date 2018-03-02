
{- 
This module is devoted to the processing of geological data.
-}

module Geometry (
  
  -- * 4D point
  Point,

  -- * 3D vector
    Vect,

  -- * geological vector
    GVect,

  -- * geological axis
    GAxis
    
 ) where


import qualified Data.Array as A
import Data.Fixed


-- Constants
p_dist_thresh = 1.0e-7
v_angle_thresh = 1.0
v_min_magn_thresh = 1.0e-9
min_vector_magn_diff = 1.0e-9


-- Utils
--------

-- | mod 360°
m360 :: Double -> Double
-- |
-- Examples:
-- >>> m360 10
-- 10.0
-- >>> m360 360
-- 0.0
-- >>> m360 (-50)
-- 310.0
-- >>> m360 400
-- 40.0
m360 v = v `mod'` 360.0

-- | convert to degrees
degrees :: Double -> Double
degrees x = x * 180.0 / pi

-- | convert to radians
radians :: Double -> Double
radians x = x * pi / 180.0

-- | calculates the angle (in degrees, positive anti-clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
angle_east_aclock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> angle_east_aclock 0 0
-- Nothing
-- >>> angle_east_aclock 1 0
-- Just 0.0
-- >>> angle_east_aclock 1 1
-- Just 45.0
-- >>> angle_east_aclock 0 5
-- Just 90.0
-- >>> angle_east_aclock (-1) 0
-- Just 180.0
-- >>> angle_east_aclock 0 (-4)
-- Just 270.0
-- >>> angle_east_aclock 1 (-1)
-- Just 315.0
angle_east_aclock x y = case (x, y) of
                          (0.0, 0.0) -> Nothing
                          (_,   _)   -> Just (m360 $ degrees $ atan2 y x)

-- | calculate the angle from East (in degrees, positive clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
angle_east_clock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> angle_east_clock 0 0
-- Nothing
-- >>> angle_east_clock 1 0
-- Just 0.0
-- >>> angle_east_clock 1 (-1)
-- Just 45.0
-- >>> angle_east_clock 0 (-4)
-- Just 90.0
-- >>> angle_east_clock (-1) 0
-- Just 180.0
-- >>> angle_east_clock 0 5
-- Just 270.0
-- >>> angle_east_clock 1 1
-- Just 315.0
angle_east_clock x y = let ang = angle_east_aclock x y
                        in case ang of
                           Nothing -> Nothing
                           Just a  -> Just (m360 $ 360.0 - a)

-- | calculate the angle from North (in degrees, positive clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
angle_north_clock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> angle_north_clock 0 0
-- Nothing
-- >>> angle_north_clock 0 5
-- Just 0.0
-- >>> angle_north_clock 1 1
-- Just 45.0
-- >>> angle_north_clock 1 0
-- Just 90.0
-- >>> angle_north_clock 1 (-1)
-- Just 135.0
-- >>> angle_north_clock 0 (-1)
-- Just 180.0
-- >>> angle_north_clock (-1) 0
-- Just 270.0
-- >>> angle_north_clock 0 (-4)
-- Just 180.0
angle_north_clock x y = let ang = angle_east_clock x y
                          in case ang of
                            Nothing -> Nothing
                            Just a  -> Just (m360 $ 90.0 + a)


-- | opposite trend
opposite_trend :: Double -> Double
-- |
-- Examples:
-- >>> opposite_trend 0
-- 180.0
-- >>> opposite_trend 45
-- 225.0
-- >>> opposite_trend 90
-- 270.0
-- >>> opposite_trend 180
-- 0.0
-- >>> opposite_trend 270
-- 90.0
opposite_trend x = m360 $ 180.0 + x


-- |  Calculates the colatitude angle from the top
plng2colatTop :: Double -> Double
-- |
-- Examples:
-- >>> plng2colatTop 90
-- 180.0
-- >>> plng2colatTop 45
-- 135.0
-- >>> plng2colatTop 0
-- 90.0
-- >>> plng2colatTop (-45)
-- 45.0
-- >>> plng2colatTop (-90)
-- 0.0
plng2colatTop plng = 90.0 + plng


-- |  Calculates the colatitude angle from the bottom
plng2colatBottom :: Double -> Double
-- |
-- Examples:
-- >>> plng2colatBottom 90
-- 0.0
-- >>> plng2colatBottom 45
-- 45.0
-- >>> plng2colatBottom 0
-- 90.0
-- >>> plng2colatBottom (-45)
-- 135.0
-- >>> plng2colatBottom (-90)
-- 180.0
plng2colatBottom plng = 90.0 - plng


-- | slope (in degrees) given horizontal and vertical lengths
-- | both input are assumed positive
slope :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> slope 0 0
-- Nothing
-- >>> slope 1 1
-- Just 45.0
-- >>> slope 1 0
-- Just 0.0
slope h v = case (h, v) of
              (0.0, 0.0) -> Nothing
              (0.0, _)   -> Just 90.0
              (_, _)     -> Just (degrees $ atan2 v h)


-- Point
--------
  
-- | Point in a 3D space
data Point = Point {px, py, pz :: Double} deriving (Show)


-- | Mapping on point
pmap :: (Double -> Double) -> Point -> Point
pmap f (Point x y z) = Point (f x) (f y) (f z)


-- | Zipping on two points
pzip :: (Double -> Double -> Double) -> Point -> Point -> Point
pzip f (Point x1 y1 z1) (Point x2 y2 z2) = Point (f x1 x2) (f y1 y2) (f z1 z2)


-- | Folding on point
pfold :: (Double -> Double -> Double) -> Point -> Double
pfold f (Point x y z) = f x (f y z)


instance Num Point where
  (+) = pzip (+)
  (-) = pzip (-)
  negate = pmap negate


-- | Conversion to array
p2a :: Point -> A.Array Int Double
p2a (Point x y z) = A.array (1, 3) [(1, x), (2, y), (3, z)]


-- | Point on the origin
pO :: Point
pO = Point 0 0 0


-- | Point on x axis
pX :: Point
pX = Point 1 0 0


-- | Point on y axis
pY :: Point
pY = Point 0 1 0


-- | Point on z axis
pZ :: Point
pZ = Point 0 0 1


-- | Projection of point on the x-y plane
p_xy :: Point -> Point
-- |
-- Examples:
-- >>> p_xy (Point 2 3 4)
-- Point {px = 2.0, py = 3.0, pz = 0.0}
p_xy (Point x y z) = Point x y 0


-- | Projection of point on the x-z plane
p_xz :: Point -> Point
-- |
-- Examples:
-- >>> p_xz (Point 2 3 4)
-- Point {px = 2.0, py = 0.0, pz = 4.0}
p_xz (Point x y z) = Point x 0 z


-- | Projection of point on the y-z plane
p_yz :: Point -> Point
-- |
-- Examples:
-- >>> p_yz (Point 2 3 4)
-- Point {px = 0.0, py = 3.0, pz = 4.0}
p_yz (Point x y z) = Point 0 y z


-- | Length (magnitude) of a point
plen :: Point -> Double
-- |
-- Examples:
-- >>> plen (Point 1 1 1)
-- 1.7320508075688772
plen (Point x y z) = sqrt(x*x + y*y + z*z)


-- | Invert point position 
pinv :: Point -> Point
-- |
-- Examples:
-- >>> pinv (Point 1 0 4)
-- Point {px = -1.0, py = -0.0, pz = -4.0} 
pinv = negate


-- | x delta between two points
pdx :: Point -> Point -> Double
pdx (Point x1 _ _) (Point x2 _ _) = x2 - x1


-- | y delta between two points
pdy :: Point -> Point -> Double
pdy (Point _ y1 _) (Point _ y2 _) = y2 - y1


-- | z delta between two points
pdz :: Point -> Point -> Double
pdz (Point _ _ z1) (Point _ _  z2) = z2 - z1


-- | 2D (horizontal) distance between two points
pd2D :: Point -> Point -> Double
pd2D (Point x1 y1 _) (Point x2 y2 _) = sqrt ((x2 - x1)^2 + (y2 - y1)^2)


-- | 3D distance between two points
pd3D :: Point -> Point -> Double
pd3D (Point x1 y1 z1) (Point x2 y2 z2) = sqrt ((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)


-- | Coincidence between two points
p_coinc :: Point -> Point -> Bool
-- |
-- Examples:
-- p_coinc (Point 1 2 3) (Point 1 2 3)
-- True
-- p_coinc (Point 1 2 3) (Point 2 2 3)
-- False
p_coinc p1 p2 = let pts_dist = pd3D p1 p2
                  in (pts_dist < p_dist_thresh)


-- | Translate point by a given amount
pshift :: Point -> Double -> Double -> Double -> Point
-- |
-- Examples:
-- pshift (Point 1 0 3) 2 7 4
-- Point {px = 3.0, py = 7.0, pz = 7.0}
pshift (Point x y z) sx sy sz = Point (x+sx) (y+sy) (z+sz)


-- | Conversion from point to vector
p2v :: Point -> Vect
-- |
-- Examples:
-- >>> p2v (Point 1 0 0)
-- Vect {x = 1.0, y = 0.0, z = 0.0}
p2v (Point x y z) = Vect x y z


-- Vector
---------

-- | Vector in a 3D cartesian space.
-- | x is East-directed, y North-directed and z upward-directed
-- | Implementation inspired to Data.Vector module
data Vect = Vect {x, y, z :: Double} deriving (Eq, Ord, Show)


-- | mapping on Vect
vmap :: (Double -> Double) -> Vect -> Vect
vmap f (Vect x y z) = Vect (f x) (f y) (f z)


-- | zipping on two Vects
vzip :: (Double -> Double -> Double) -> Vect -> Vect -> Vect
vzip f (Vect x1 y1 z1) (Vect x2 y2 z2) = Vect (f x1 x2) (f y1 y2) (f z1 z2)


-- | folding on Vect
vfold :: (Double -> Double -> Double) -> Vect -> Double
vfold f (Vect x y z) = f x (f y z)


-- | vector dot product
vdot :: Vect -> Vect -> Double
-- |
-- Examples:
-- >>> vdot (Vect 1 2 4) (Vect 0 3 15)
-- 66.0
-- >>> vdot (Vect 2 0 3) (Vect 3 2 7)
-- 27.0
-- >>> vdot (Vect 1 0 0) (Vect 1 0 0)
-- 1.0
-- >>> vdot (Vect 1 0 0) (Vect 0 1 0)
-- 0.0
-- >>> vdot (Vect 1 0 0) (Vect (-1) 0 0)
-- -1.0
vdot v1 v2 = vfold (+) $ vzip (*) v1 v2


-- | vector cross product
vcross :: Vect -> Vect -> Vect
-- |
-- Examples:
-- >>> vcross (Vect 1 0 0) (Vect 0 1 0)
-- Vect {x = 0.0, y = 0.0, z = 1.0}
-- >>> vcross vectZ vectX
-- Vect {x = 0.0, y = 1.0, z = 0.0}
vcross (Vect x1 y1 z1) (Vect x2 y2 z2) = Vect {
  x = y1 * z2 - y2 * z1,
  y = z1 * x2 - z2 * x1,
  z = x1 * y2 - x2 * y1}


instance Num Vect where
  (+) = vzip (+)
  (-) = vzip (-)
  (*) = vcross
  negate = vmap negate


-- | conversion to array
v2a :: Vect -> A.Array Int Double
v2a (Vect x y z) = A.array (1, 3) [(1, x), (2, y), (3, z)]


-- | versor parallel to x axis
vectX :: Vect
vectX = Vect 1 0 0


-- | versor parallel to y axis
vectY :: Vect
vectY = Vect 0 1 0


-- | versor parallel to z axis
vectZ :: Vect
vectZ = Vect 0 0 1


-- | horizontal length of a vector
vlenh :: Vect -> Double
-- |
-- Examples:
-- >>> vlenh (Vect 0 1 23)
-- 1.0
-- >>> vlenh (Vect 1 1 1)
-- 1.4142135623730951
-- >>> vlenh (Vect 3 4 2)
-- 5.0
vlenh (Vect x y _) = sqrt (x*x + y*y)


-- | length (magnitude) of a vector
vlen :: Vect -> Double
-- |
-- Examples:
-- >>> vlen (Vect 1 1 1)
-- 1.7320508075688772
vlen (Vect x y z) = sqrt(x*x + y*y + z*z)


-- | is vector with almost zero components
v_almost_zero :: Vect -> Bool
-- |
-- Examples:
-- >>> v_almost_zero (Vect 1.0e-10 1.0e-10 1.0e-10)
-- True
-- >>> v_almost_zero (Vect 1 0 2)
-- False
v_almost_zero v = let l = vlen v
                   in (l < v_min_magn_thresh)


-- | is regular vector
v_notzero :: Vect -> Bool
-- |
-- Examples:
-- >>> v_notzero (Vect 1 0 0)
-- True
-- >>> v_notzero (Vect 0 0 0)
-- False
v_notzero v = not (v_almost_zero v)


-- | is vector almost unitary
v_almost_unit :: Vect -> Bool
-- |
-- Examples:
-- >>> v_almost_unit (Vect 1 0 0)
-- True
-- >>> v_almost_unit (Vect 0 2 7)
-- False
v_almost_unit v = abs (1 - (vlen v)) < min_vector_magn_diff


-- | are two vectors almost equal
v_almost_equal :: Vect -> Vect -> Bool
-- |
-- Examples:
-- >>> v_almost_equal (Vect 1 2 0) (Vect 1 2 1.0e-12)
-- True
-- >>> v_almost_equal (Vect 0.1 3.4 0.7) (Vect 6.2 5.3 9.2)
-- False
v_almost_equal v1 v2 = let l1 = vlen v1
                           l2 = vlen v2
                        in
                         (abs (l1 - l2)) < min_vector_magn_diff
                         

-- | invert vector direction 
vinv :: Vect -> Vect
-- |
-- Examples:
-- >>> vinv (Vect 1 0 4)
-- Vect {x = -1.0, y = -0.0, z = -4.0} 
vinv = negate


-- | vector multiplication by scalar
vmul :: Vect -> Double -> Vect
-- |
-- Examples:
-- >>> vmul (Vect 1 0 3) 4
-- Vect {x = 4.0, y = 0.0, z = 12.0}
-- >>> vmul (Vect 1 3 2.5) 0.5
-- Vect {x = 0.5, y = 1.5, z = 1.25}
vmul v s = vmap (*s) v


-- | vector division by scalar
(//) :: Vect -> Double -> Maybe Vect
-- |
-- Examples:
-- >>> (Vect 3 7 11) // 2
-- Just (Vect {x = 1.5, y = 3.5, z = 5.5})
-- >>> (Vect 3 4 9) // 0.0
-- Nothing
(//) v s = case s of
             0.0 -> Nothing
             _   -> Just (vmap (/s) v)


-- | vector normalization
versor :: Vect -> Maybe Vect
-- |
-- Examples:
-- >>> versor (Vect 2 0 0)
-- Just (Vect {x = 1.0, y = 0.0, z = 0.0}) 
-- >>> versor (Vect 0 0 0)
-- Nothing
versor v = case l of
             0.0 -> Nothing
             _  -> v // l
           where l = vlen v


-- | check if upward-pointing vector
v_is_upward :: Vect -> Bool
-- |
-- Examples:
-- >>> v_is_upward (Vect 1 0 3) 
-- True
-- >>> v_is_upward (Vect 0 (-2) (-3))
-- False
v_is_upward (Vect _ _ z) = z > 0.0


-- | check if downward-pointing vector
v_is_downward :: Vect -> Bool
-- |
-- Examples:
-- >>> v_is_downward (Vect 1 0 (-3)) 
-- True
-- >>> v_is_downward (Vect 0 (-2) 3)
-- False
v_is_downward (Vect _ _ z) = z < 0.0


-- | upward-pointing vector
vup :: Vect -> Vect
-- |
-- Examples:
-- >>> vup (Vect 1 1 1)
-- Vect {x = 1.0, y = 1.0, z = 1.0}
-- >>> vup (Vect (-1) (-1) (-1))
-- Vect {x = 1.0, y = 1.0, z = 1.0}
vup (Vect x y z) =
  if (z < 0.0)
    then vinv (Vect x y z)
    else Vect x y z


-- | calculate a new downward-pointing vector
vdown :: Vect -> Vect
-- |
-- Examples:
-- >>> vdown (Vect 1 1 1)
-- Vect {x = -1.0, y = -1.0, z = -1.0}
-- >>> vdown (Vect 3 (-7) (-1))
-- Vect {x = 3.0, y = -7.0, z = -1.0}
vdown (Vect x y z) =
  if (z > 0.0)
    then vinv (Vect x y z)
    else Vect x y z


-- | trend of a vector
--   (degrees, clockwise from North, range 0°-360°)
vtrend :: Vect -> Maybe Double
-- |
-- Examples:
-- >>> vtrend (Vect 1 0 0)
-- Just 90.0
-- >>> vtrend (Vect 0 1 0)
-- Just 0.0
-- >>> vtrend (Vect 1 1 0)
-- Just 45.0
-- >>> vtrend (Vect 1 (-1) 0)
-- Just 135.0 
-- >>> vtrend (Vect 0 (-1) 0)
-- Just 180.0
-- >>> vtrend (Vect (-1) (-1) 0)
-- Just 225.0
-- >>> vtrend (Vect (-1) 0 0)
-- Just 270.0
-- >>> vtrend (Vect (-1) 1 0)
-- Just 315.0
-- >>> vtrend (Vect 1 1 10)
-- Just 45.0
vtrend (Vect x y z) = angle_north_clock x y


-- | slope of vector
--   (degrees, positive: downward-directed, negative: upward-dir., range -90°/90°
vslope :: Vect -> Maybe Double
-- |
-- Examples:
-- >>> vslope (Vect 1 0 (-1))
-- Just 45.0
-- >>> vslope (Vect 1 0 1)
-- Just (-45.0)
-- >>> vslope (Vect 0 1 0)
-- Just 0.0
-- >>> vslope (Vect 0 0 1)
-- Just (-90.0)
-- >>> vslope (Vect 0 0 (-1))
-- Just 90.0
-- >>> vslope (Vect 0 0 0)
-- Nothing
vslope v = let h = vlenh v
               zv = z v
               sl = slope h (abs zv)
            in case sl of
              Nothing -> Nothing
              Just slp   -> if (zv <= 0.0)
                             then Just slp
                             else Just (-slp)  

              
-- | calculate the geological vector parallel to the Vect instance.
-- | Trend range: [0°, 360°[
-- | Plunge range: [-90°, 90°], with negative values for upward-pointing
-- | geological axes and positive values for downward-pointing axes.
vToGVect :: Vect -> Maybe GVect
-- |
-- Examples:
-- >>> vToGVect (Vect 1 1 1)
-- Just (GVect {tr = 45.0, pl = -35.264389682754654})
-- >>> vToGVect (Vect 0 1 1)
-- Just (GVect {tr = 0.0, pl = -45.0})
-- >>> vToGVect (Vect 1 0 1)
-- Just (GVect {tr = 90.0, pl = -45.0})
-- >>> vToGVect (Vect 0 0 1)
-- Just (GVect {tr = 0.0, pl = -90.0})
-- >>> vToGVect (Vect 0 0 (-1))
-- Just (GVect {tr = 0.0, pl = 90.0})
-- >>> vToGVect (Vect (-1) 0 0)
-- Just (GVect {tr = 270.0, pl = 0.0})
-- >>> vToGVect (Vect 0 (-1) 0)
-- Just (GVect {tr = 180.0, pl = 0.0})
-- >>> vToGVect (Vect (-1) (-1) 0)
-- Just (GVect {tr = 225.0, pl = 0.0})
vToGVect v = let trend = vtrend v
                 plunge = vslope v
            in case (trend, plunge) of
                (_,       Just    90.0)  -> Just (GVect 0.0   90.0 )
                (_,       Just (-90.0))  -> Just (GVect 0.0 (-90.0))
                (_,       Nothing     )  -> Nothing
                (Nothing,    _        )  -> Nothing               
                (Just t,  Just p      )  -> Just (GVect t p)


-- | calculate the geological axis parallel to the Vect instance.
-- | Trend range: [0°, 360°[
-- | Plunge range: [-90°, 90°], with negative values for upward-pointing
-- | geological axes and positive values for downward-pointing axes.
gaxis :: Vect -> Maybe GAxis
-- |
-- | Examples:
-- | >>> gaxis (Vect 0 1 1)
-- | GAxis(000.00, -45.00)
-- | >>> gaxis (Vect 1 0 1)
-- | GAxis(090.00, -45.00)
-- | >>> gaxis (Vect 0 0 1)
-- | GAxis(000.00, -90.00)
-- | >>> gaxis (Vect 0 0 (-1))
-- | GAxis(000.00, +90.00)
-- | >>> gaxis (Vect (-1) 0, 0)
-- | GAxis(270.00, +00.00)
-- | >>> gaxis (Vect 0 (-1) 0)
-- | GAxis(180.00, +00.00)
-- | >>> gaxis (Vect (-1) (-1) 0)
-- | GAxis(225.00, +00.00)
-- | >>> gaxis (Vect 0 0 0)
-- | Nothing
gaxis v = gvToGa $ vToGVect v


-- | Cosine of the angle between two versors.
versors_cos_angle :: Vect -> Vect -> Double
-- |
-- Examples:
-- >>> versors_cos_angle (Vect 1 0 0) (Vect 0 0 1)
-- 0.0
-- >>> versors_cos_angle (Vect 1 0 0) (Vect (-1) 0 0)
-- -1.0
-- >>> versors_cos_angle (Vect 1 0 0) (Vect 1 0 0)
-- 1.0
versors_cos_angle v1 v2 = let  dp = vdot v1 v2 in
                            if (dp > 1.0)
                              then 1.0
                            else if (dp < (-1.0))
                              then -1.0
                            else dp
                                             

-- | Return the cosine of the angle between two vectors.
cos_angle :: Vect -> Vect -> Maybe Double
-- |
-- Examples:
-- >>> cos_angle (Vect 1 0 0) (Vect 0 0 1)
-- Just 0.0
-- >>> cos_angle (Vect 1 0 0) (Vect (-1) 0 0)
-- Just (-1.0)
-- >>> cos_angle (Vect 1 0 0) (Vect 1 0 0)
-- Just 1.0
-- >>> cos_angle (Vect 0 0 0) (Vect 1 0 0)
-- Nothing
-- >>> cos_angle (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
cos_angle v1 v2 = let uv1 = versor v1
                      uv2 = versor v2
                   in case (uv1, uv2) of
                     (Nothing, _)         -> Nothing
                     (_, Nothing)         -> Nothing
                     (Just vv1, Just vv2) -> Just (versors_cos_angle vv1 vv2)


-- | Calculate angle between two vectors, as degrees
-- | in 0° - 180° range.
vAngle :: Vect -> Vect -> Maybe Double
-- |
-- Examples:
-- >>> vAngle (Vect 1 0 0) (Vect 0 0 1)
-- Just 90.0
-- >>> vAngle (Vect 1 0 0) (Vect (-1) 0 0)
-- Just 180.0
-- >>> vAngle (Vect 0 0 1) (Vect 0 0 (-1))
-- Just 180.0
-- >>> vAngle (Vect 1 1 1) (Vect 1 1 1)
-- Just 0.0
-- >>> vAngle (Vect 0 0 0) (Vect 1 0 0)
-- Nothing
-- >>> vAngle (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
vAngle v1 v2 = let cos_ang = cos_angle v1 v2
               in case cos_ang of
                 Nothing -> Nothing
                 Just ca -> Just (degrees $ acos ca)


-- | Determine whether two vectors are sub-parallel
v_subparallel :: Vect -> Vect -> Maybe Bool
-- |  
-- Examples:
-- >>> v_subparallel (Vect 1 0 0) (Vect 1 0 0)
-- Just True
-- >>> v_subparallel (Vect 1 0 0) (Vect 0 0 1)
-- Just False
-- >>> v_subparallel (Vect 1 0 0) (Vect (-1) 0 0)
-- Just False
-- >>> v_subparallel (Vect 0 0 0) (Vect 1 0 0)
-- Nothing
-- >>> v_subparallel (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
v_subparallel v1 v2 = let ang = vAngle v1 v2
                       in case ang of
                         Nothing  -> Nothing
                         Just a   -> Just (a < v_angle_thresh)

                                     
-- | Check whether two vectors are sub-orhogonal
v_is_suborthogonal :: Vect -> Vect -> Maybe Bool
-- |
-- Example:
-- >>> v_is_suborthogonal (Vect 1 0 0) (Vect 1 1 0)
-- Just False
-- >>> v_is_suborthogonal (Vect 1 0 0) (Vect 0 1 0)
-- Just True
-- >>> v_is_suborthogonal (Vect 1 0 0) (Vect 0 1 1)
-- Just True
-- >>> v_is_suborthogonal (Vect 1 0 0) (Vect 0 0.9999999999999 0)
-- Just True
-- >>> v_is_suborthogonal (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
v_is_suborthogonal v1 v2 = let ang = vAngle v1 v2
                            in case ang of
                              Nothing -> Nothing
                              Just a  -> Just ((90.0 - a) < v_angle_thresh)

-- | Matrix multiplication of a vector
-- TODO



{- | GVect

Geological vector.
Defined by trend and plunge (both in degrees):
 - trend: [0.0, 360.0[ clockwise, from 0 (North):
 - plunge: [-90.0, 90.0] - negative value: upward pointing axis,
                           positive values: downward axis;
-}

data GVect = GVect {tr, pl :: Double} deriving (Eq, Ord, Show)


-- | mapping on GVect
gvmap :: (Double -> Double) -> GVect -> GVect
gvmap f (GVect tr pl) = GVect (f tr) (f pl)


-- | zipping on GVect
gvzip :: (Double -> Double -> Double) -> GVect -> GVect -> GVect
gvzip f (GVect tr1 pl1) (GVect tr2 pl2) = GVect (f tr1 tr2) (f pl1 pl2)


-- | Return trend of the geological direction
-- | Unit is degree, range is [0, 360[
gvTr :: GVect -> Double
-- |
-- Example:
-- >>> gvTr (GVect 420 (-17))
-- 60.0
-- >>> gvTr (GVect (-20) 49)
-- 340.0
gvTr (GVect tr _) = m360 tr


-- | Return plunge of the geological direction
-- | Unit is degree, range is [-90, 90]
-- | Negative values are upward pointing
-- | Positive values are downward pointing
gvPl :: GVect -> Double
-- |
-- Example:
-- >>> gvPl (GVect 420 (-17))
-- -17.0
gvPl (GVect _ pl) = pl


-- | Return trend and plunge of the geological direction
gvTp :: GVect -> (Double, Double)
-- |
-- Example:
-- >>> gvTp (GVect (-90) (-45))
-- (270.0,-45.0)
gvTp (GVect tr pl) = (m360 tr, pl)


-- | Calculates the colatitude from the North (top)
-- | return an angle bewtween 0 and 180 (as degrees)
gvColatitudeNorth :: GVect -> Double
-- |
-- Examples:
-- >>> gvColatitudeNorth (GVect 320 90)
-- 180.0
-- >>> gvColatitudeNorth (GVect 320 45)
-- 135.0
-- >>> gvColatitudeNorth (GVect 320 0)
-- 90.0
-- >>> gvColatitudeNorth (GVect 320 (-45))
-- 45.0
-- >>> gvColatitudeNorth (GVect 320 (-90))
-- 0.0
gvColatitudeNorth (GVect _ pl) = plng2colatTop pl


-- | Calculates the colatitude from the South (bottom)
-- | return an angle bewtween 0 and 180 (as degrees)
gvColatitudeSouth :: GVect -> Double
-- |
-- Examples:
-- >>> gvColatitudeSouth (GVect 320 90)
-- 0.0
-- >>> gvColatitudeSouth (GVect 320 45)
-- 45.0
-- >>> gvColatitudeSouth (GVect 320 0)
-- 90.0
-- >>> gvColatitudeSouth (GVect 320 (-45))
-- 135.0
-- >>> gvColatitudeSouth (GVect 320 (-90))
-- 180.0
gvColatitudeSouth (GVect _ pl) = plng2colatBottom pl


-- | Return a copy of the GVect instance
gvCopy :: GVect -> GVect
-- |
-- Example:
-- >>> gvCopy (GVect 10 20)
-- GVect {tr = 10.0, pl = 20.0}
gvCopy (GVect tr pl) = GVect tr pl


-- | GVect opposite
gvOpposite :: GVect -> GVect
-- |
-- Example:
-- >>> gvOpposite (GVect 30 20)
-- GVect {tr = 210.0, pl = -20.0}
-- >>> gvOpposite (GVect 190 (-40))
-- GVect {tr = 10.0, pl = 40.0}
gvOpposite (GVect tr0 pl0) =
  let tr1 = opposite_trend tr0
      pl1 = - pl0
   in GVect tr1 pl1


-- | Determine the GVect instance that
-- | is normal and downward-sliping
gvNormalDown :: GVect -> GVect
-- |
-- Examples:
-- >>> gvNormalDown (GVect 90 0)
-- GVect {tr = 270.0, pl = 90.0}
-- >>> gvNormalDown (GVect 45 30)
-- GVect {tr = 225.0, pl = 60.0}
-- >>> gvNormalDown (GVect 180 90)
-- GVect {tr = 0.0, pl = 0.0}
gvNormalDown (GVect tr pl) =
  let dip_dir = opposite_trend tr
      dip_angle = 90.0 - pl
   in GVect dip_dir dip_angle


-- | Unit vector from GVect instance
gvToVersor :: GVect -> Vect
gvToVersor (GVect trend plunge) =
  let tr = radians trend
      pl = radians plunge
      cos_tr = cos tr
      sin_tr = sin tr
      cos_pl = cos pl
      sin_pl = sin pl
      north = cos_pl * cos_tr
      east  = cos_pl * sin_tr
      down  = sin_pl
   in Vect east north (-down)

-- | Conversion to GPlane
gvToGPlane :: GVect -> GPlane
-- |

-- Examples:
-- >>> gvToGPlane (GVect 90 45)
-- GPlane {az = 90.0, dip = 45.0}
gvToGPlane (GVect tr pl) = GPlane tr pl

  
-- |  Check whether the instance is pointing upward
gvIsUpward :: GVect -> Bool
-- |
--    Examples:
--    >>> gvIsUpward (GVect 10 15)
--    False
--    >>> gvIsUpward (GVect 257.4 0.0)
--    False
--    >>> gvIsUpward (GVect 90 (-45))
--    True
gvIsUpward gv = (z $ gvToVersor gv) > 0.0


-- |  Check whether the instance is pointing downward
gvIsDownward :: GVect -> Bool
-- |
-- Examples:
-- >>> gvIsDownward (GVect 10 15)
-- True
-- >>> gvIsDownward (GVect 257.4 0.0)
-- False
-- >>> gvIsDownward (GVect 90 (-45))
-- False
gvIsDownward gv = (z $ gvToVersor gv) < 0.0


-- | Return upward-point geological vector
gvUp :: GVect -> GVect
-- |
-- Examples:
-- >>> gvUp (GVect 90 (-45))
-- GVect {tr = 90.0, pl = -45.0}
-- >>> gvUp (GVect 180 45)
-- GVect {tr = 0.0, pl = -45.0}
-- >>> gvUp (GVect 0 0)
-- GVect {tr = 0.0, pl = 0.0}
-- >>> gvUp (GVect 0 90)
-- GVect {tr = 180.0, pl = -90.0}
-- >>> gvUp (GVect 0 (-90))
-- GVect {tr = 0.0, pl = -90.0}
gvUp gv = if (gvIsDownward gv)
           then gvOpposite gv
           else gvCopy gv
                

-- | Return downward-point geological vector
gvDown :: GVect -> GVect
-- |
-- Examples:
-- >>> gvDown (GVect 90 (-45))
-- GVect {tr = 270.0, pl = 45.0}
-- >>> gvDown (GVect 180 45)
-- GVect {tr = 180.0, pl = 45.0}
-- >>> gvDown (GVect 0 0)
-- GVect {tr = 0.0, pl = 0.0}
-- >>> gvDown (GVect 0 90)
-- GVect {tr = 0.0, pl = 90.0}
-- >>> gvDown (GVect 0 (-90))
-- GVect {tr = 180.0, pl = 90.0}
gvDown gv = if (gvIsUpward gv)
           then gvOpposite gv
           else gvCopy gv


-- |  Calculate angle (in degrees) between the two GVect instances.
-- |  Range is 0°-180°.
gvAngle :: GVect -> GVect -> Double
-- |
-- Examples:
-- >>> gvAngle (GVect 0 90) (GVect 90 0)
-- 90.0
-- >>> gvAngle (GVect 0 0) (GVect 270 0)
-- 90.00000000000001
-- >>> gvAngle (GVect 0 0) (GVect 0 0)
-- 0.0
-- >>> gvAngle (GVect 0 0) (GVect 180 0)
-- 180.0
-- >>> gvAngle (GVect 90 0) (GVect 270 0)
-- 180.0
gvAngle gv1 gv2 = degrees $ acos $ versors_cos_angle (gvToVersor gv1) (gvToVersor gv2)


-- |  Check that two GVect are sub-parallel
gvAlmostParallel :: GVect -> GVect -> Bool
-- |
-- Examples:
-- >>> gvAlmostParallel (GVect 0 90) (GVect 90 0)
-- False
-- >>> gvAlmostParallel (GVect 0 0) (GVect 0 1e-6)
-- True
-- >>> gvAlmostParallel (GVect 0 90) (GVect 180 0)
-- False
-- >>> gvAlmostParallel (GVect 0 90) (GVect 0 (-90))
-- False
gvAlmostParallel gv1 gv2 = (gvAngle gv1 gv2) < v_angle_thresh


-- | Return the geological plane normal to the geological vector.
gvNormalGPlane :: GVect -> GPlane
-- |
-- Examples:
-- >>> gvNormalGPlane (GVect 0 45)
-- GPlane {az = 180.0, dip = 45.0}
-- >>> gvNormalGPlane (GVect 0 (-45))
-- GPlane {az = 0.0, dip = 45.0}
-- >>> gvNormalGPlane (GVect 0 90)
-- GPlane {az = 180.0, dip = 0.0}
gvNormalGPlane gv =
  let down_axis = gvDown gv
      norm_down = gvNormalDown down_axis
   in gvToGPlane norm_down
       

-- | Plane common to two GVect instances
gvCommonPlane :: GVect -> GVect -> Maybe GPlane
-- |
-- Examples:
-- >>> gvCommonPlane (GVect 0 0) (GVect 90 0)
-- Just (GPlane {az = 180.0, dip = 0.0})
-- >>> gvCommonPlane (GVect 0 0) (GVect 90 90)
-- Just (GPlane {az = 90.0, dip = 90.0})
-- >>> gvCommonPlane (GVect 45 0) (GVect 135 45)
-- Just (GPlane {az = 135.0, dip = 44.99999999999999})
-- >>> gvCommonPlane (GVect 315 45) (GVect 135 45)
-- Just (GPlane {az = 225.0, dip = 90.0})
-- >>> gvCommonPlane (GVect 0 0) (GVect 90 0)
-- Just (GPlane {az = 180.0, dip = 0.0})
-- >>> gvCommonPlane (GVect 0 0) (GVect 90 90)
-- Just (GPlane {az = 90.0, dip = 90.0})
-- >>> gvCommonPlane (GVect 45 0) (GVect 135 45)
-- Just (GPlane {az = 135.0, dip = 44.99999999999999})
-- >>> gvCommonPlane (GVect 315 45) (GVect 135 45)
-- Just (GPlane {az = 225.0, dip = 90.0})
-- >>> gvCommonPlane (GVect 0 0) (GVect 0 0)
-- Nothing
gvCommonPlane gv1 gv2 = let v1 = gvToVersor gv1
                            v2 = gvToVersor gv2
                            normal_v = vcross v1 v2
                            normal_gv = vToGVect normal_v
                         in case normal_gv of
                           Nothing -> Nothing
                           Just gv -> Just (gvNormalGPlane gv)


{- Geological axis. While GAxis is non-directional, the geological vector (GVect) is directional.
    Defined by trend and plunge (both in degrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
-}

data GAxis = GAxis {t, p :: Double} deriving (Eq, Ord, Show)

gvToGa gv = case gv of
             Nothing -> Nothing
             Just (GVect tr pl) -> Just (GAxis tr pl)


{- Geological plane.
    Defined by dip direction azimuth and dip angle (both in degrees):
     - dip direction azimuth: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
-}

data GPlane = GPlane {az, dip :: Double} deriving (Eq, Ord, Show)




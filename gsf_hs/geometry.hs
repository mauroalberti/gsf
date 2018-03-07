
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
kPDistThresh = 1.0e-7
kVAngleThresh = 1.0
kVMinMagnThresh = 1.0e-9
kVMinMagnDiff = 1.0e-9


-- Utils
--------

-- | mod 360°
uM360 :: Double -> Double
-- |
-- Examples:
-- >>> uM360 10
-- 10.0
-- >>> uM360 360
-- 0.0
-- >>> uM360 (-50)
-- 310.0
-- >>> uM360 400
-- 40.0
uM360 v = v `mod'` 360.0

-- | convert to uDegrees
uDegrees :: Double -> Double
uDegrees x = x * 180.0 / pi

-- | convert to uRadians
uRadians :: Double -> Double
uRadians x = x * pi / 180.0

-- | calculates the angle (in uDegrees, positive anti-clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
uAngleEastAclock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> uAngleEastAclock 0 0
-- Nothing
-- >>> uAngleEastAclock 1 0
-- Just 0.0
-- >>> uAngleEastAclock 1 1
-- Just 45.0
-- >>> uAngleEastAclock 0 5
-- Just 90.0
-- >>> uAngleEastAclock (-1) 0
-- Just 180.0
-- >>> uAngleEastAclock 0 (-4)
-- Just 270.0
-- >>> uAngleEastAclock 1 (-1)
-- Just 315.0
uAngleEastAclock x y = case (x, y) of
                          (0.0, 0.0) -> Nothing
                          (_,   _)   -> Just (uM360 $ uDegrees $ atan2 y x)

-- | calculate the angle from East (in uDegrees, positive clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
uAngleEastClock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> uAngleEastClock 0 0
-- Nothing
-- >>> uAngleEastClock 1 0
-- Just 0.0
-- >>> uAngleEastClock 1 (-1)
-- Just 45.0
-- >>> uAngleEastClock 0 (-4)
-- Just 90.0
-- >>> uAngleEastClock (-1) 0
-- Just 180.0
-- >>> uAngleEastClock 0 5
-- Just 270.0
-- >>> uAngleEastClock 1 1
-- Just 315.0
uAngleEastClock x y = let ang = uAngleEastAclock x y
                        in case ang of
                           Nothing -> Nothing
                           Just a  -> Just (uM360 $ 360.0 - a)

-- | calculate the angle from North (in uDegrees, positive clockwise) in the 2D plane
-- given x (East) and y (North) Cartesian components
uAngleNorthClock :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> uAngleNorthClock 0 0
-- Nothing
-- >>> uAngleNorthClock 0 5
-- Just 0.0
-- >>> uAngleNorthClock 1 1
-- Just 45.0
-- >>> uAngleNorthClock 1 0
-- Just 90.0
-- >>> uAngleNorthClock 1 (-1)
-- Just 135.0
-- >>> uAngleNorthClock 0 (-1)
-- Just 180.0
-- >>> uAngleNorthClock (-1) 0
-- Just 270.0
-- >>> uAngleNorthClock 0 (-4)
-- Just 180.0
uAngleNorthClock x y = let ang = uAngleEastClock x y
                          in case ang of
                            Nothing -> Nothing
                            Just a  -> Just (uM360 $ 90.0 + a)


-- | opposite trend
uOppositeTrend :: Double -> Double
-- |
-- Examples:
-- >>> uOppositeTrend 0
-- 180.0
-- >>> uOppositeTrend 45
-- 225.0
-- >>> uOppositeTrend 90
-- 270.0
-- >>> uOppositeTrend 180
-- 0.0
-- >>> uOppositeTrend 270
-- 90.0
uOppositeTrend x = uM360 $ 180.0 + x


-- |  Calculates the colatitude angle from the top
uPlngToColatTop :: Double -> Double
-- |
-- Examples:
-- >>> uPlngToColatTop 90
-- 180.0
-- >>> uPlngToColatTop 45
-- 135.0
-- >>> uPlngToColatTop 0
-- 90.0
-- >>> uPlngToColatTop (-45)
-- 45.0
-- >>> uPlngToColatTop (-90)
-- 0.0
uPlngToColatTop plng = 90.0 + plng


-- |  Calculates the colatitude angle from the bottom
uPlngToColatBottom :: Double -> Double
-- |
-- Examples:
-- >>> uPlngToColatBottom 90
-- 0.0
-- >>> uPlngToColatBottom 45
-- 45.0
-- >>> uPlngToColatBottom 0
-- 90.0
-- >>> uPlngToColatBottom (-45)
-- 135.0
-- >>> uPlngToColatBottom (-90)
-- 180.0
uPlngToColatBottom plng = 90.0 - plng


-- | uSlope (in uDegrees) given horizontal and vertical lengths
-- | both input are assumed positive
uSlope :: Double -> Double -> Maybe Double
-- |
-- Examples:
-- >>> uSlope 0 0
-- Nothing
-- >>> uSlope 1 1
-- Just 45.0
-- >>> uSlope 1 0
-- Just 0.0
uSlope h v = case (h, v) of
              (0.0, 0.0) -> Nothing
              (0.0, _)   -> Just 90.0
              (_, _)     -> Just (uDegrees $ atan2 v h)


-- Point
--------
  
-- | Point in a 3D space
data Point = Point {px, py, pz :: Double} deriving (Show)


-- | Mapping on point
pMap :: (Double -> Double) -> Point -> Point
pMap f (Point x y z) = Point (f x) (f y) (f z)


-- | Zipping on two points
pZip :: (Double -> Double -> Double) -> Point -> Point -> Point
pZip f (Point x1 y1 z1) (Point x2 y2 z2) = Point (f x1 x2) (f y1 y2) (f z1 z2)


-- | Folding on point
pFold :: (Double -> Double -> Double) -> Point -> Double
pFold f (Point x y z) = f x (f y z)


instance Num Point where
  (+) = pZip (+)
  (-) = pZip (-)
  negate = pMap negate


-- | Conversion to array
pToArray :: Point -> A.Array Int Double
pToArray (Point x y z) = A.array (1, 3) [(1, x), (2, y), (3, z)]


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
pXY :: Point -> Point
-- |
-- Examples:
-- >>> pXY (Point 2 3 4)
-- Point {px = 2.0, py = 3.0, pz = 0.0}
pXY (Point x y z) = Point x y 0


-- | Projection of point on the x-z plane
pXZ :: Point -> Point
-- |
-- Examples:
-- >>> pXZ (Point 2 3 4)
-- Point {px = 2.0, py = 0.0, pz = 4.0}
pXZ (Point x y z) = Point x 0 z


-- | Projection of point on the y-z plane
pYZ :: Point -> Point
-- |
-- Examples:
-- >>> pYZ (Point 2 3 4)
-- Point {px = 0.0, py = 3.0, pz = 4.0}
pYZ (Point x y z) = Point 0 y z


-- | Length (magnitude) of a point
pLen :: Point -> Double
-- |
-- Examples:
-- >>> pLen (Point 1 1 1)
-- 1.7320508075688772
pLen (Point x y z) = sqrt(x*x + y*y + z*z)


-- | Invert point position 
pInvert :: Point -> Point
-- |
-- Examples:
-- >>> pInvert (Point 1 0 4)
-- Point {px = -1.0, py = -0.0, pz = -4.0} 
pInvert = negate


-- | x delta between two points
pDeltaX :: Point -> Point -> Double
pDeltaX (Point x1 _ _) (Point x2 _ _) = x2 - x1


-- | y delta between two points
pDeltaY :: Point -> Point -> Double
pDeltaY (Point _ y1 _) (Point _ y2 _) = y2 - y1


-- | z delta between two points
pDeltaZ :: Point -> Point -> Double
pDeltaZ (Point _ _ z1) (Point _ _  z2) = z2 - z1


-- | 2D (horizontal) distance between two points
pDist2D :: Point -> Point -> Double
pDist2D (Point x1 y1 _) (Point x2 y2 _) = sqrt ((x2 - x1)^2 + (y2 - y1)^2)


-- | 3D distance between two points
pDist3D :: Point -> Point -> Double
pDist3D (Point x1 y1 z1) (Point x2 y2 z2) = sqrt ((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)


-- | Coincidence between two points
pAreCoinc :: Point -> Point -> Bool
-- |
-- Examples:
-- pAreCoinc (Point 1 2 3) (Point 1 2 3)
-- True
-- pAreCoinc (Point 1 2 3) (Point 2 2 3)
-- False
pAreCoinc p1 p2 = let pts_dist = pDist3D p1 p2
                  in (pts_dist < kPDistThresh)


-- | Translate point by a given amount
pShift :: Point -> Double -> Double -> Double -> Point
-- |
-- Examples:
-- pShift (Point 1 0 3) 2 7 4
-- Point {px = 3.0, py = 7.0, pz = 7.0}
pShift (Point x y z) sx sy sz = Point (x+sx) (y+sy) (z+sz)


-- | Conversion from point to vector
pToVect :: Point -> Vect
-- |
-- Examples:
-- >>> pToVect (Point 1 0 0)
-- Vect {x = 1.0, y = 0.0, z = 0.0}
pToVect (Point x y z) = Vect x y z


-- Vector
---------

-- | Vector in a 3D cartesian space.
-- | x is East-directed, y North-directed and z upward-directed
-- | Implementation inspired to Data.Vector module
data Vect = Vect {x, y, z :: Double} deriving (Eq, Ord, Show)


-- | mapping on Vect
vMap :: (Double -> Double) -> Vect -> Vect
vMap f (Vect x y z) = Vect (f x) (f y) (f z)


-- | zipping on two Vects
vZip :: (Double -> Double -> Double) -> Vect -> Vect -> Vect
vZip f (Vect x1 y1 z1) (Vect x2 y2 z2) = Vect (f x1 x2) (f y1 y2) (f z1 z2)


-- | folding on Vect
vFold :: (Double -> Double -> Double) -> Vect -> Double
vFold f (Vect x y z) = f x (f y z)


-- | vector dot product
vDot :: Vect -> Vect -> Double
-- |
-- Examples:
-- >>> vDot (Vect 1 2 4) (Vect 0 3 15)
-- 66.0
-- >>> vDot (Vect 2 0 3) (Vect 3 2 7)
-- 27.0
-- >>> vDot (Vect 1 0 0) (Vect 1 0 0)
-- 1.0
-- >>> vDot (Vect 1 0 0) (Vect 0 1 0)
-- 0.0
-- >>> vDot (Vect 1 0 0) (Vect (-1) 0 0)
-- -1.0
vDot v1 v2 = vFold (+) $ vZip (*) v1 v2


-- | vector cross product
vCross :: Vect -> Vect -> Vect
-- |
-- Examples:
-- >>> vCross (Vect 1 0 0) (Vect 0 1 0)
-- Vect {x = 0.0, y = 0.0, z = 1.0}
-- >>> vCross vectZ vectX
-- Vect {x = 0.0, y = 1.0, z = 0.0}
vCross (Vect x1 y1 z1) (Vect x2 y2 z2) = Vect {
  x = y1 * z2 - y2 * z1,
  y = z1 * x2 - z2 * x1,
  z = x1 * y2 - x2 * y1}


instance Num Vect where
  (+) = vZip (+)
  (-) = vZip (-)
  (*) = vCross
  negate = vMap negate


-- | conversion to array
vToArray :: Vect -> A.Array Int Double
vToArray (Vect x y z) = A.array (1, 3) [(1, x), (2, y), (3, z)]


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
vLenH :: Vect -> Double
-- |
-- Examples:
-- >>> vLenH (Vect 0 1 23)
-- 1.0
-- >>> vLenH (Vect 1 1 1)
-- 1.4142135623730951
-- >>> vLenH (Vect 3 4 2)
-- 5.0
vLenH (Vect x y _) = sqrt (x*x + y*y)


-- | length (magnitude) of a vector
vLen :: Vect -> Double
-- |
-- Examples:
-- >>> vLen (Vect 1 1 1)
-- 1.7320508075688772
vLen (Vect x y z) = sqrt(x*x + y*y + z*z)


-- | is vector with almost zero components
vIsAlmostZero :: Vect -> Bool
-- |
-- Examples:
-- >>> vIsAlmostZero (Vect 1.0e-10 1.0e-10 1.0e-10)
-- True
-- >>> vIsAlmostZero (Vect 1 0 2)
-- False
vIsAlmostZero v = let l = vLen v
                   in (l < kVMinMagnThresh)


-- | is regular vector
vIsNotZero :: Vect -> Bool
-- |
-- Examples:
-- >>> vIsNotZero (Vect 1 0 0)
-- True
-- >>> vIsNotZero (Vect 0 0 0)
-- False
vIsNotZero v = not (vIsAlmostZero v)


-- | is vector almost unitary
vIsAlmostUnit :: Vect -> Bool
-- |
-- Examples:
-- >>> vIsAlmostUnit (Vect 1 0 0)
-- True
-- >>> vIsAlmostUnit (Vect 0 2 7)
-- False
vIsAlmostUnit v = abs (1 - (vLen v)) < kVMinMagnDiff


-- | are two vectors almost equal
vAreAlmostEqual :: Vect -> Vect -> Bool
-- |
-- Examples:
-- >>> vAreAlmostEqual (Vect 1 2 0) (Vect 1 2 1.0e-12)
-- True
-- >>> vAreAlmostEqual (Vect 0.1 3.4 0.7) (Vect 6.2 5.3 9.2)
-- False
vAreAlmostEqual v1 v2 = let l1 = vLen v1
                            l2 = vLen v2
                        in
                         (abs (l1 - l2)) < kVMinMagnDiff
                         

-- | invert vector direction 
vInvert :: Vect -> Vect
-- |
-- Examples:
-- >>> vInvert (Vect 1 0 4)
-- Vect {x = -1.0, y = -0.0, z = -4.0} 
vInvert = negate


-- | vector multiplication by scalar
vMul :: Vect -> Double -> Vect
-- |
-- Examples:
-- >>> vMul (Vect 1 0 3) 4
-- Vect {x = 4.0, y = 0.0, z = 12.0}
-- >>> vMul (Vect 1 3 2.5) 0.5
-- Vect {x = 0.5, y = 1.5, z = 1.25}
vMul v s = vMap (*s) v


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
             _   -> Just (vMap (/s) v)


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
           where l = vLen v


-- | check if upward-pointing vector
vIsUpward :: Vect -> Bool
-- |
-- Examples:
-- >>> vIsUpward (Vect 1 0 3) 
-- True
-- >>> vIsUpward (Vect 0 (-2) (-3))
-- False
vIsUpward (Vect _ _ z) = z > 0.0


-- | check if downward-pointing vector
vIsDownward :: Vect -> Bool
-- |
-- Examples:
-- >>> vIsDownward (Vect 1 0 (-3)) 
-- True
-- >>> vIsDownward (Vect 0 (-2) 3)
-- False
vIsDownward (Vect _ _ z) = z < 0.0


-- | upward-pointing vector
vUpward :: Vect -> Vect
-- |
-- Examples:
-- >>> vUpward (Vect 1 1 1)
-- Vect {x = 1.0, y = 1.0, z = 1.0}
-- >>> vUpward (Vect (-1) (-1) (-1))
-- Vect {x = 1.0, y = 1.0, z = 1.0}
vUpward (Vect x y z) =
  if (z < 0.0)
    then vInvert (Vect x y z)
    else Vect x y z


-- | calculate a new downward-pointing vector
vDownward :: Vect -> Vect
-- |
-- Examples:
-- >>> vDownward (Vect 1 1 1)
-- Vect {x = -1.0, y = -1.0, z = -1.0}
-- >>> vDownward (Vect 3 (-7) (-1))
-- Vect {x = 3.0, y = -7.0, z = -1.0}
vDownward (Vect x y z) =
  if (z > 0.0)
    then vInvert (Vect x y z)
    else Vect x y z


-- | trend of a vector
--   (uDegrees, clockwise from North, range 0°-360°)
vTrend :: Vect -> Maybe Double
-- |
-- Examples:
-- >>> vTrend (Vect 1 0 0)
-- Just 90.0
-- >>> vTrend (Vect 0 1 0)
-- Just 0.0
-- >>> vTrend (Vect 1 1 0)
-- Just 45.0
-- >>> vTrend (Vect 1 (-1) 0)
-- Just 135.0 
-- >>> vTrend (Vect 0 (-1) 0)
-- Just 180.0
-- >>> vTrend (Vect (-1) (-1) 0)
-- Just 225.0
-- >>> vTrend (Vect (-1) 0 0)
-- Just 270.0
-- >>> vTrend (Vect (-1) 1 0)
-- Just 315.0
-- >>> vTrend (Vect 1 1 10)
-- Just 45.0
vTrend (Vect x y z) = uAngleNorthClock x y


-- | uSlope of vector
--   (uDegrees, positive: downward-directed, negative: upward-dir., range -90°/90°
vSlope :: Vect -> Maybe Double
-- |
-- Examples:
-- >>> vSlope (Vect 1 0 (-1))
-- Just 45.0
-- >>> vSlope (Vect 1 0 1)
-- Just (-45.0)
-- >>> vSlope (Vect 0 1 0)
-- Just 0.0
-- >>> vSlope (Vect 0 0 1)
-- Just (-90.0)
-- >>> vSlope (Vect 0 0 (-1))
-- Just 90.0
-- >>> vSlope (Vect 0 0 0)
-- Nothing
vSlope v = let h = vLenH v
               zv = z v
               sl = uSlope h (abs zv)
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
vToGVect v = let trend = vTrend v
                 plunge = vSlope v
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
vToGAxis :: Vect -> Maybe GAxis
-- |
-- | Examples:
-- | >>> vToGAxis (Vect 0 1 1)
-- | GAxis(000.00, -45.00)
-- | >>> vToGAxis (Vect 1 0 1)
-- | GAxis(090.00, -45.00)
-- | >>> vToGAxis (Vect 0 0 1)
-- | GAxis(000.00, -90.00)
-- | >>> vToGAxis (Vect 0 0 (-1))
-- | GAxis(000.00, +90.00)
-- | >>> vToGAxis (Vect (-1) 0, 0)
-- | GAxis(270.00, +00.00)
-- | >>> vToGAxis (Vect 0 (-1) 0)
-- | GAxis(180.00, +00.00)
-- | >>> vToGAxis (Vect (-1) (-1) 0)
-- | GAxis(225.00, +00.00)
-- | >>> vToGAxis (Vect 0 0 0)
-- | Nothing
vToGAxis v = let gv = vToGVect v
              in case gv of
               Nothing -> Nothing
               Just val_gv -> Just (gvToGAxis val_gv)


-- | Cosine of the angle between two versors.
versCosAngle :: Vect -> Vect -> Double
-- |
-- Examples:
-- >>> versCosAngle (Vect 1 0 0) (Vect 0 0 1)
-- 0.0
-- >>> versCosAngle (Vect 1 0 0) (Vect (-1) 0 0)
-- -1.0
-- >>> versCosAngle (Vect 1 0 0) (Vect 1 0 0)
-- 1.0
versCosAngle v1 v2 = let  dp = vDot v1 v2 in
                            if (dp > 1.0)
                              then 1.0
                            else if (dp < (-1.0))
                              then -1.0
                            else dp
                                             

-- | Return the cosine of the angle between two vectors.
vCosAngle :: Vect -> Vect -> Maybe Double
-- |
-- Examples:
-- >>> vCosAngle (Vect 1 0 0) (Vect 0 0 1)
-- Just 0.0
-- >>> vCosAngle (Vect 1 0 0) (Vect (-1) 0 0)
-- Just (-1.0)
-- >>> vCosAngle (Vect 1 0 0) (Vect 1 0 0)
-- Just 1.0
-- >>> vCosAngle (Vect 0 0 0) (Vect 1 0 0)
-- Nothing
-- >>> vCosAngle (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
vCosAngle v1 v2 = let uv1 = versor v1
                      uv2 = versor v2
                   in case (uv1, uv2) of
                     (Nothing, _)         -> Nothing
                     (_, Nothing)         -> Nothing
                     (Just vv1, Just vv2) -> Just (versCosAngle vv1 vv2)


-- | Calculate angle between two vectors, as uDegrees
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
vAngle v1 v2 = let cos_ang = vCosAngle v1 v2
               in case cos_ang of
                 Nothing -> Nothing
                 Just ca -> Just (uDegrees $ acos ca)


-- | Determine whether two vectors are sub-parallel
vAreSubParallel :: Vect -> Vect -> Maybe Bool
-- |  
-- Examples:
-- >>> vAreSubParallel (Vect 1 0 0) (Vect 1 0 0)
-- Just True
-- >>> vAreSubParallel (Vect 1 0 0) (Vect 0 0 1)
-- Just False
-- >>> vAreSubParallel (Vect 1 0 0) (Vect (-1) 0 0)
-- Just False
-- >>> vAreSubParallel (Vect 0 0 0) (Vect 1 0 0)
-- Nothing
-- >>> vAreSubParallel (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
vAreSubParallel v1 v2 = let ang = vAngle v1 v2
                       in case ang of
                         Nothing  -> Nothing
                         Just a   -> Just (a < kVAngleThresh)

                                     
-- | Check whether two vectors are sub-orhogonal
vAreSubOrthogonal :: Vect -> Vect -> Maybe Bool
-- |
-- Example:
-- >>> vAreSubOrthogonal (Vect 1 0 0) (Vect 1 1 0)
-- Just False
-- >>> vAreSubOrthogonal (Vect 1 0 0) (Vect 0 1 0)
-- Just True
-- >>> vAreSubOrthogonal (Vect 1 0 0) (Vect 0 1 1)
-- Just True
-- >>> vAreSubOrthogonal (Vect 1 0 0) (Vect 0 0.9999999999999 0)
-- Just True
-- >>> vAreSubOrthogonal (Vect 1 0 0) (Vect 0 0 0)
-- Nothing
vAreSubOrthogonal v1 v2 = let ang = vAngle v1 v2
                            in case ang of
                              Nothing -> Nothing
                              Just a  -> Just ((90.0 - a) < kVAngleThresh)

-- | Matrix multiplication of a vector
-- TODO



{- | GVect

Geological vector.
Defined by trend and plunge (both in uDegrees):
 - trend: [0.0, 360.0[ clockwise, from 0 (North):
 - plunge: [-90.0, 90.0] - negative value: upward pointing axis,
                           positive values: downward axis;
-}

data GVect = GVect {tr, pl :: Double} deriving (Eq, Ord, Show)


-- | mapping on GVect
gvMap :: (Double -> Double) -> GVect -> GVect
gvMap f (GVect tr pl) = GVect (f tr) (f pl)


-- | zipping on GVect
gvZip :: (Double -> Double -> Double) -> GVect -> GVect -> GVect
gvZip f (GVect tr1 pl1) (GVect tr2 pl2) = GVect (f tr1 tr2) (f pl1 pl2)


-- | Return trend of the geological direction
-- | Unit is degree, range is [0, 360[
gvTr :: GVect -> Double
-- |
-- Example:
-- >>> gvTr (GVect 420 (-17))
-- 60.0
-- >>> gvTr (GVect (-20) 49)
-- 340.0
gvTr (GVect tr _) = uM360 tr


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
gvTrPl :: GVect -> (Double, Double)
-- |
-- Example:
-- >>> gvTrPl (GVect (-90) (-45))
-- (270.0,-45.0)
gvTrPl (GVect tr pl) = (uM360 tr, pl)


-- | Calculates the colatitude from the North (top)
-- | return an angle bewtween 0 and 180 (as uDegrees)
gvColatNorth :: GVect -> Double
-- |
-- Examples:
-- >>> gvColatNorth (GVect 320 90)
-- 180.0
-- >>> gvColatNorth (GVect 320 45)
-- 135.0
-- >>> gvColatNorth (GVect 320 0)
-- 90.0
-- >>> gvColatNorth (GVect 320 (-45))
-- 45.0
-- >>> gvColatNorth (GVect 320 (-90))
-- 0.0
gvColatNorth (GVect _ pl) = uPlngToColatTop pl


-- | Calculates the colatitude from the South (bottom)
-- | return an angle bewtween 0 and 180 (as uDegrees)
gvColatSouth :: GVect -> Double
-- |
-- Examples:
-- >>> gvColatSouth (GVect 320 90)
-- 0.0
-- >>> gvColatSouth (GVect 320 45)
-- 45.0
-- >>> gvColatSouth (GVect 320 0)
-- 90.0
-- >>> gvColatSouth (GVect 320 (-45))
-- 135.0
-- >>> gvColatSouth (GVect 320 (-90))
-- 180.0
gvColatSouth (GVect _ pl) = uPlngToColatBottom pl


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
  let tr1 = uOppositeTrend tr0
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
  let dip_dir = uOppositeTrend tr
      dip_angle = 90.0 - pl
   in GVect dip_dir dip_angle


-- | Unit vector from GVect instance
gvToVersor :: GVect -> Vect
gvToVersor (GVect trend plunge) =
  let tr = uRadians trend
      pl = uRadians plunge
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
gvUpward :: GVect -> GVect
-- |
-- Examples:
-- >>> gvUpward (GVect 90 (-45))
-- GVect {tr = 90.0, pl = -45.0}
-- >>> gvUpward (GVect 180 45)
-- GVect {tr = 0.0, pl = -45.0}
-- >>> gvUpward (GVect 0 0)
-- GVect {tr = 0.0, pl = 0.0}
-- >>> gvUpward (GVect 0 90)
-- GVect {tr = 180.0, pl = -90.0}
-- >>> gvUpward (GVect 0 (-90))
-- GVect {tr = 0.0, pl = -90.0}
gvUpward gv = if (gvIsDownward gv)
           then gvOpposite gv
           else gvCopy gv
                

-- | Return downward-point geological vector
gvDownward :: GVect -> GVect
-- |
-- Examples:
-- >>> gvDownward (GVect 90 (-45))
-- GVect {tr = 270.0, pl = 45.0}
-- >>> gvDownward (GVect 180 45)
-- GVect {tr = 180.0, pl = 45.0}
-- >>> gvDownward (GVect 0 0)
-- GVect {tr = 0.0, pl = 0.0}
-- >>> gvDownward (GVect 0 90)
-- GVect {tr = 0.0, pl = 90.0}
-- >>> gvDownward (GVect 0 (-90))
-- GVect {tr = 180.0, pl = 90.0}
gvDownward gv = if (gvIsUpward gv)
           then gvOpposite gv
           else gvCopy gv


-- |  Calculate angle (in uDegrees) between the two GVect instances.
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
gvAngle gv1 gv2 = uDegrees $ acos $ versCosAngle (gvToVersor gv1) (gvToVersor gv2)


-- |  Check that two GVect are sub-parallel
gvAreAlmostParallel :: GVect -> GVect -> Bool
-- |
-- Examples:
-- >>> gvAreAlmostParallel (GVect 0 90) (GVect 90 0)
-- False
-- >>> gvAreAlmostParallel (GVect 0 0) (GVect 0 1e-6)
-- True
-- >>> gvAreAlmostParallel (GVect 0 90) (GVect 180 0)
-- False
-- >>> gvAreAlmostParallel (GVect 0 90) (GVect 0 (-90))
-- False
gvAreAlmostParallel gv1 gv2 = (gvAngle gv1 gv2) < kVAngleThresh


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
  let down_axis = gvDownward gv
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
                            normal_v = vCross v1 v2
                            normal_gv = vToGVect normal_v
                         in case normal_gv of
                           Nothing -> Nothing
                           Just gv -> Just (gvNormalGPlane gv)

-- | Converts a geological vector to a geological axis
gvToGAxis :: GVect -> GAxis
-- |
-- >>> gvToGAxis (GVect 220 32)
-- GAxis {t = 220.0, p = 32.0}
gvToGAxis (GVect tr pl) = GAxis tr pl


{- Calculate the GVect instance that is normal to the two provided sources.
Angle between sources must be larger than kVAngleThresh,
otherwise Nothing will be returned.
-}
gvNormalGVect :: GVect -> GVect -> Maybe GVect
-- |
-- Examples:
-- >>> gvNormalGVect (GVect 0 0) (GVect 0 0.5)
-- Nothing
-- >>> gvNormalGVect (GVect 0 0) (GVect 179.1 0)
-- Nothing
-- >>> gvNormalGVect (GVect 0 0) (GVect 5.1 0)
-- Just (GVect {tr = 0.0, pl = 90.0})
-- >>> gvNormalGVect (GVect 90 45) (GVect 90 0)
-- Just (GVect {tr = 180.0, pl = 0.0})
gvNormalGVect gv1 gv2 = let angle = gvAngle gv1 gv2
                            delta_angle = if (angle < kVAngleThresh) then False
                                          else if (angle >= (180.0 - kVAngleThresh)) then False
                                          else True
                            in case delta_angle of
                             False -> Nothing
                             True  -> vToGVect $ vCross (gvToVersor gv1) (gvToVersor gv2)


{- Geological axis. While GAxis is non-directional, the geological vector (GVect) is directional.
    Defined by trend and plunge (both in uDegrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
-}

data GAxis = GAxis {t, p :: Double} deriving (Eq, Ord, Show)

-- | Converts a geological axis to a geological vector
gaToGVect :: GAxis -> GVect
-- |
-- Example:
-- >>> gaToGVect (GAxis 220 32)
-- GVect {tr = 220.0, pl = 32.0}
gaToGVect (GAxis t p) = GVect t p


-- | Converts a geological axis to a versor
gaToVers :: GAxis -> Vect
-- |
-- Examples:
-- >>> gaToVers (GAxis 90 0)
-- Vect {x = 1.0, y = 6.123233995736766e-17, z = -0.0}
-- >>> gaToVers (GAxis 0 45)
-- Vect {x = 0.0, y = 0.7071067811865476, z = -0.7071067811865475}
-- >>> gaToVers (GAxis 0 90)
-- Vect {x = 0.0, y = 6.123233995736766e-17, z = -1.0}
-- >>> gaToVers (GAxis 270 (-90))
-- Vect {x = -6.123233995736766e-17, y = -1.1248198369963932e-32, z = 1.0}
gaToVers ga = gvToVersor $ gaToGVect ga


-- | Calculate angle (in uDegrees) between the two GAxis instances.
-- | Range: 0.0 - 90.0
gaAngle :: GAxis -> GAxis -> Double
-- |
-- Examples:
-- >>> gaAngle (GAxis 0 90) (GAxis 90 0)
-- 90.0
-- >>> gaAngle (GAxis 0 0) (GAxis 270 0)
-- 89.99999999999999
-- >>> gaAngle (GAxis 0 0) (GAxis 0 0)
-- 0.0
-- >>> gaAngle (GAxis 0 0) (GAxis 180 0)
-- 0.0
-- >>> gaAngle (GAxis 0 0) (GAxis 179 0)
-- 0.9999999999998863
-- >>> gaAngle (GAxis 0 (-90)) (GAxis 0 90)
-- 0.0
-- >>> gaAngle (GAxis 90 0) (GAxis 315 0)
-- 44.99999999999997
gaAngle ga1 ga2 = let angle_vers = uDegrees $ acos $ versCosAngle (gaToVers ga1) (gaToVers ga2)
                   in min  angle_vers (180.0 - angle_vers)

-- | Check that two GAxis are sub-parallel
gaAreSubParallel :: GAxis -> GAxis -> Bool
-- |
-- Examples:
-- >>> gaAreSubParallel (GAxis 0 90) (GAxis 90 0)
-- False
-- >>> gaAreSubParallel (GAxis 0 0) (GAxis 0 1.0e-6)
-- True
-- >>> gaAreSubParallel (GAxis 0 0) (GAxis 180 0)
-- True
-- >>> gaAreSubParallel (GAxis 90 0) (GAxis 270 0)
-- True
-- >>> gaAreSubParallel (GAxis 0 90) (GAxis 0 (-90))
-- True
gaAreSubParallel ga1 ga2 = (gaAngle ga1 ga2) <= kVAngleThresh
  

-- | Check that two GAxis are sub-orthogonal
gaAreSubOrthogonal :: GAxis -> GAxis -> Bool
-- |
-- Examples:
-- >>> gaAreSubOrthogonal (GAxis 0 90) (GAxis 90 0)
-- True
-- >>> gaAreSubOrthogonal (GAxis 0 0) (GAxis 0 1.0e-6)
-- False
-- >>> gaAreSubOrthogonal (GAxis 0 0) (GAxis 180 0)
-- False
-- >>> gaAreSubOrthogonal (GAxis 90 0) (GAxis 270 89.5)
-- True
-- >>> gaAreSubOrthogonal (GAxis 0 90) (GAxis 0 0.5)
-- True
gaAreSubOrthogonal ga1 ga2 = (90.0 - (gaAngle ga1 ga2)) < kVAngleThresh


-- | Calculate the geological plane normal to a given GAxis instance.
gaNormalGPlane :: GAxis -> GPlane
-- |
-- Example:
-- >>> gaNormalGPlane (GAxis 0 90)
-- GPlane {az = 180.0, dip = 0.0}
-- >>> gaNormalGPlane (GAxis 0 0)
-- GPlane {az = 180.0, dip = 90.0}
-- >>> gaNormalGPlane (GAxis 45 45)
-- GPlane {az = 225.0, dip = 45.0}
gaNormalGPlane ga = gvNormalGPlane $ gaToGVect ga


{- Geological plane.
    Defined by dip direction azimuth and dip angle (both in uDegrees):
     - dip direction azimuth: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
-}

data GPlane = GPlane {az, dip :: Double} deriving (Eq, Ord, Show)




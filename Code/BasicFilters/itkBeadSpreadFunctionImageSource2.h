/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef _itkBeadSpreadFunctionImageSource2_h
#define _itkBeadSpreadFunctionImageSource2_h

#include "itkCommand.h"
#include "itkEllipseSpatialObject.h"
#include "itkParametricImageSource.h"
#include "itkPointSet.h"
#include "itkPointSetConvolutionImageFilter.h"
#include "itkShiftScaleImageFilter.h"

namespace itk
{

/** \class BeadSpreadFunctionImageSource2
 *
 * \brief Generates a synthetic bead-spread function that is the
 * convolution of a sphere with a ParametricImageSource that generates
 * a convolution kernel. This version uses the
 * PointSetConvolutionImageFilter on a set of points generated by
 * rasterizing a sphere.
 *
 * \ingroup DataSources Multithreaded
*/
template < class TOutputImage >
class ITK_EXPORT BeadSpreadFunctionImageSource2 :
    public ParametricImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef BeadSpreadFunctionImageSource2        Self;
  typedef ParametricImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Typedef for output types. */
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename PointType::VectorType           VectorType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::IndexValueType IndexValueType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef ParametricImageSource< TOutputImage >
    KernelImageSourceType;
  typedef typename KernelImageSourceType::Pointer
    KernelImageSourcePointer;
  typedef PointSet< PixelType, TOutputImage::ImageDimension >
    PointSetType;
  typedef PointSetConvolutionImageFilter< PointSetType, TOutputImage, TOutputImage >
    ConvolverType;
  typedef typename ConvolverType::Pointer
    ConvolverPointer;
  typedef ShiftScaleImageFilter< TOutputImage, TOutputImage >
    RescaleImageFilterType;
  typedef typename RescaleImageFilterType::Pointer
    RescaleImageFilterPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BeadSpreadFunctionImageSource2,ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Set/get the size of the output image. */
  void SetSize(const SizeType & size);
  itkGetConstReferenceMacro(Size, SizeType);

  /** Set/get the spacing of the output image (in nanometers). */
  void SetSpacing(const SpacingType & spacing);
  const SpacingType & GetSpacing() const;

  /** Set/get the origin of the output image (in nanometers). */
  virtual void SetOrigin(const PointType & origin);
  const PointType & GetOrigin() const;

  /** Set/get the bead radius (in nanometers). */
  //itkSetMacro(BeadRadius, double);
  void SetBeadRadius( double radius );
  itkGetConstMacro(BeadRadius, double);

  /** Set/get the shear in the X direction. */
  void SetShearX(double shear);
  double GetShearX() const;

  /** Set/get the shear in the Y direction. */
  void SetShearY(double shear);
  double GetShearY() const;

  /** Set/get the bead center. */
  itkSetMacro(BeadCenter, PointType);
  itkGetConstReferenceMacro(BeadCenter, PointType);

  /** Set/get the background value. */
  itkSetMacro(IntensityShift, double);
  itkGetConstMacro(IntensityShift, double);

  /** Set/get the maximum intensity. */
  itkSetMacro(IntensityScale, double);
  itkGetConstMacro(IntensityScale, double);

  /** Set/get the convolution kernel source. */
  virtual void SetKernelSource( KernelImageSourceType* source );
  itkGetObjectMacro(KernelSource, KernelImageSourceType);

  /** Set/get a single parameter value. */
  virtual void SetParameter(unsigned int index, ParametersValueType value);
  virtual ParametersValueType GetParameter(unsigned int index) const;

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Gets the number of bead-spread function parameters. */
  virtual unsigned int GetNumberOfBeadSpreadFunctionParameters() const;

  /** Get/set the z-coordinate of the image z-plane at the given index. */
  void SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int);

  /** Get/set use of custom z coordinates. */
  void SetUseCustomZCoordinates(bool use);
  bool GetUseCustomZCoordinates();

  /** Callback evoked whenever the KernelSource is modified. */
  virtual void KernelModified();

protected:
  BeadSpreadFunctionImageSource2();
  virtual ~BeadSpreadFunctionImageSource2();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Recompute the PointSet resulting from rasterizing the bead
   * shape. */
  void RasterizeShape();

  /** This class is implicitly multi-threaded because its member filters
   * are mulithreaded, so we go with a "single-threaded"
   * implementation here. */
  virtual void GenerateData();

  virtual void GenerateOutputInformation();

private:
  BeadSpreadFunctionImageSource2(const BeadSpreadFunctionImageSource2&); // purposely not implemented
  void operator=(const BeadSpreadFunctionImageSource2&); // purposely not implemented

  /** This is needed because the GetSize() method in the
   * PointSetConvolutionImageFilter returns a copy of the Size, not a
   * reference. Because this class derives from ParametricImageSource,
   * it must override its GetSize() method to return a
   * reference. Returning the copy from the
   * PointSetConvolutionImageFilter isn't safe because it returns a
   * reference to a temporary variable allocated on the stack. */
  SizeType m_Size;

  double      m_BeadRadius;
  PointType   m_BeadCenter;
  SpacingType m_BeadSampleSpacing;

  double m_IntensityShift; // Additive background constant
  double m_IntensityScale; // The maximum intensity value

  KernelImageSourcePointer   m_KernelSource;
  ConvolverPointer           m_Convolver;
  RescaleImageFilterPointer  m_RescaleFilter;

  typedef SimpleMemberCommand< Self > MemberCommandType;
  typedef typename MemberCommandType::Pointer MemberCommandPointer;

  MemberCommandPointer m_ModifiedEventCommand;
  unsigned long        m_ObserverTag;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBeadSpreadFunctionImageSource2.txx"
#endif


#endif // _itkBeadSpreadFunctionImageSource2_h

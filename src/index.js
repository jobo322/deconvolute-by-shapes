'use strict';

// const { signalToPeaks } = require('nmr-processing');
const { getShape1D } = require('ml-peak-shape-generator');

const { DefaultParameters } = require('./DefaultParameters');
const { assert } = require('./assert');
//I will optimise signals (clusters of peaks) if the J

const properties = ['init', 'min', 'max', 'gradientDifference'];

const signalModel = {
  x: 0,
  y: 1,
  shape: { kind: 'pseudoVoigt' },
  parameters: { x: { init: 0, max: 1, min: -1, gradientDifference: 0.1 } },
  pattern: [{ x: -1, y: 1 }, { x: 1, y: 1 }] //A doublet
}

const internalSignalModel = {

}

function getSumOfShapes(internalSignals) {
  return function sumOfShapes(parameters) {
    return (currentX) => {
      let totalY = 0;
      for (const signal of internalSignals) {
        const delta = parameters[signal.fromIndex];
        const intensity = parameters[signal.fromIndex + 1];
        for (let i = 2; i <= signal.toIndex; i++) {
          //@ts-expect-error Not simply to solve the issue
          signal.shapeFct[signal.parameters[i]] = parameters[signal.fromIndex + i];
        }
        for (let peak of signal.pattern) {
          const { x, y } = peak;
          totalY += y * intensity * signal.shapeFct.fct(currentX - x - delta);
        }
      }
      return totalY;
    };
  };
}

export function getInternalPeaks(
  peaks,
  minMaxY,
  options,
) {
  let index = 0;
  let internalPeaks = [];
  for (const peak of peaks) {
    const { pattern = { x: 0, y: 1 } } = peak;
    const shape = peak.shape
      ? peak.shape
      : options.shape
        ? options.shape
        : { kind: 'gaussian' };

    const shapeFct = getShape1D(shape);

    //@ts-expect-error Should disappear with next release of peak-shape-generator
    const parameters = ['x', 'y', ...shapeFct.getParameters()];

    const propertiesValues = {
      min: [],
      max: [],
      init: [],
      gradientDifference: [],
    };

    for (let parameter of parameters) {
      for (let property of properties) {
        // check if the property is specified in the peak
        let propertyValue = getParameterByKey(
          parameter,
          property,
          peak,
        );
        if (propertyValue) {
          propertyValue = getNormalizedValue(
            propertyValue,
            parameter,
            property,
            minMaxY,
          );

          propertiesValues[property].push(propertyValue);
          continue;
        }
        // check if there are some global option, it could be a number or a callback

        let generalParameterValue = getParameterByKey(
          parameter,
          property,
          options,
        );
        if (generalParameterValue) {
          if (typeof generalParameterValue === 'number') {
            generalParameterValue = getNormalizedValue(
              generalParameterValue,
              parameter,
              property,
              minMaxY,
            );
            propertiesValues[property].push(generalParameterValue);
            continue;
          } else {
            let value = generalParameterValue(peak);
            value = getNormalizedValue(value, parameter, property, minMaxY);
            propertiesValues[property].push(value);
            continue;
          }
        }

        // we just need to take the default parameters
        assert(
          DefaultParameters[parameter],
          `No default parameter for ${parameter}`,
        );
        const defaultParameterValues = DefaultParameters[parameter][property];
        //@ts-expect-error should never happen
        propertiesValues[property].push(defaultParameterValues(peak, shapeFct));
      }
    }

    const fromIndex = index;
    const toIndex = fromIndex + parameters.length - 1;
    index += toIndex - fromIndex + 1;

    internalPeaks.push({
      shape,
      pattern,
      shapeFct,
      parameters,
      propertiesValues,
      fromIndex,
      toIndex,
    });
  }
  return internalPeaks;
}

function getNormalizedValue(
  value,
  parameter,
  property,
  minMaxY,
) {
  if (parameter === 'y') {
    if (property === 'gradientDifference') {
      return value / minMaxY.range;
    } else {
      return (value - minMaxY.min) / minMaxY.range;
    }
  }
  return value;
}

function getParameterByKey(
  parameterKey,
  property,
  options,
) {
  return options?.parameters?.[parameterKey]?.[property];
}